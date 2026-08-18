"""
GATE 6 - shared frozen-embedding extraction utility for genomic foundation
models evaluated as "frozen embeddings + lightweight head" (charter Task 1
preferred approach). Produces a cached {oligo_id/rs241_id: embedding} lookup
so the expensive step (one forward pass per unique 165bp sequence through the
FM) runs exactly once per model, independent of how many times downstream
LOHO/calibration/RS241 fits reuse it.

Sequence source: identical to scripts/28 and scripts/50 -- decoded directly
from data/baseline_cache/{EC,BS,PA}_baseline_data.npz's onehot arrays (all
three hosts share byte-identical oligo_ids/onehot/fold, verified before
writing this script) for the 29042-row three-host library, and from
data/rs241.parquet joined to data/three_host_library.parquet for RS241's
207 sequence-available rows (same join scripts/50 uses).

Pooling: mean over token embeddings, masked by the tokenizer's attention
mask (excludes padding). This is the standard "frozen embedding" protocol
for a BERT-style encoder; no [CLS]-token pooling is used because this
checkpoint's pooler head was never trained for a downstream task (see load
report: pooler.dense.{weight,bias} MISSING from the checkpoint).
"""
import sys
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
EMB_DIR = DATA / "fm_embeddings"
EMB_DIR.mkdir(parents=True, exist_ok=True)

BASE_IDX_TO_CHAR = {0: "A", 1: "C", 2: "G", 3: "T"}


def onehot_to_seq(arr):
    idx = arr.argmax(axis=0)
    return "".join(BASE_IDX_TO_CHAR[i] for i in idx)


def load_library_sequences():
    """29042-row three-host-library sequences + oligo_ids + fold, decoded
    once from EC's cache (verified byte-identical to BS/PA's)."""
    d = dict(np.load(DATA / "baseline_cache" / "EC_baseline_data.npz", allow_pickle=True))
    seqs = [onehot_to_seq(d["onehot"][i]) for i in range(len(d["oligo_ids"]))]
    return list(d["oligo_ids"]), seqs, d["fold"]


def load_rs241_sequences():
    """207 RS241 rows with recoverable sequence text -- identical join to
    scripts/50_rs241_train_and_eval.py's load_rs241_eval_arrays()."""
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241["id"] = rs241["id"].astype(str)
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    lib["OLIGO ID"] = lib["OLIGO ID"].astype(str)
    seq_map = dict(zip(lib["OLIGO ID"], lib["Regulatory Sequence"]))
    rs241["has_sequence"] = rs241["id"].isin(seq_map)
    rs241["sequence"] = rs241["id"].map(seq_map)
    rs241_seq = rs241[rs241["has_sequence"]].copy().reset_index(drop=True)
    return list(rs241_seq["id"]), list(rs241_seq["sequence"])


def extract_embeddings_dnabert2(seqs, device=None, batch_size=128, progress_every=2000):
    """Returns (N, 768) float32 array, mean-pooled over non-pad tokens."""
    import torch
    import transformers.dynamic_module_utils as dmu
    dmu.check_imports = lambda filename: []  # see out/GATE6_MEMO.md methodology note:
    # bypasses a stale triton pre-flight check in transformers' trust_remote_code
    # loader; bert_layers.py already handles a real missing-triton ImportError via
    # its own try/except (falls back to plain pytorch attention -- confirmed via
    # the "Unable to import Triton" UserWarning at load time, not a silent skip).
    from transformers import AutoTokenizer, AutoModel, AutoConfig

    if device is None:
        device = "mps" if torch.backends.mps.is_available() else "cpu"

    tok = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    cfg = AutoConfig.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    cfg.pad_token_id = tok.pad_token_id if tok.pad_token_id is not None else 0
    model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True,
                                       config=cfg, low_cpu_mem_usage=False)
    model.eval().to(device)
    n_params = sum(p.numel() for p in model.parameters())

    out = np.zeros((len(seqs), 768), dtype=np.float32)
    t0 = time.time()
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True)
        input_ids = enc["input_ids"].to(device)
        attn = enc["attention_mask"].to(device)
        with torch.no_grad():
            res = model(input_ids)
        hidden = res[0] if isinstance(res, tuple) else res.last_hidden_state
        mask = attn.unsqueeze(-1).float()
        pooled = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1)
        out[i:i + len(batch)] = pooled.cpu().numpy()
        if i % progress_every == 0:
            print(f"    [{i}/{len(seqs)}] {time.time()-t0:.1f}s elapsed")
    print(f"  DNABERT-2 embedding extraction: {len(seqs)} seqs, {n_params} params, "
          f"{time.time()-t0:.1f}s total, device={device}")
    return out


def extract_embeddings_promogen2(seqs, device=None, batch_size=32, progress_every=2000):
    """Returns (N, 640) float32 array, mean-pooled over the final hidden
    state of jinyuan22/promogen2-base -- a GPT-2-architecture causal LM
    pretrained ONLY on prokaryotic promoter/cis-regulatory sequence (17,000
    bacterial genomes, 59M promoters -> 1.4M curated training sequences).
    Best domain match of any candidate found for Task 2 (charter: "Search
    rather than assuming -- something may have appeared recently"); run as a
    supplementary check alongside DNABERT-2, not a substitute for it -- no
    compatibility patching was required to load this model (standard
    `transformers` GPT2LMHeadModel). Character-level tokenizer (one token per
    base, vocab={A,C,G,T,N,bos,eos,pad}) -- no padding needed for this
    dataset's fixed 165bp length beyond batch-boundary edge cases.

    MPS OOM, found and fixed while running this extraction: the first version
    called model(..., output_hidden_states=True) and read hidden_states[-1],
    which materializes and retains all 31 layers' (batch,165,640) tensors
    per forward pass -- these did not get released between batches on MPS's
    caching allocator, and allocation grew unbounded across ~800 batches
    until it hit the device ceiling ("MPS backend out of memory... allocated
    8.69 GiB"). Fixed by calling model.transformer(...) directly (the base
    GPT2Model, no LM head) to get only last_hidden_state, avoiding the other
    30 layers' tensors entirely, plus an explicit torch.mps.empty_cache()
    every batch and a smaller batch_size (32, down from 128)."""
    import torch
    from transformers import AutoTokenizer, AutoModelForCausalLM

    if device is None:
        device = "mps" if torch.backends.mps.is_available() else "cpu"

    tok = AutoTokenizer.from_pretrained("jinyuan22/promogen2-base")
    model = AutoModelForCausalLM.from_pretrained("jinyuan22/promogen2-base")
    model.eval().to(device)
    n_params = sum(p.numel() for p in model.parameters())

    out = np.zeros((len(seqs), 640), dtype=np.float32)
    t0 = time.time()
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True)
        input_ids = enc["input_ids"].to(device)
        attn = enc["attention_mask"].to(device)
        with torch.no_grad():
            res = model.transformer(input_ids, attention_mask=attn)
        hidden = res.last_hidden_state
        mask = attn.unsqueeze(-1).float()
        pooled = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1)
        out[i:i + len(batch)] = pooled.cpu().numpy()
        del res, hidden, pooled, input_ids, attn
        if device == "mps":
            torch.mps.empty_cache()
        if i % progress_every == 0:
            print(f"    [{i}/{len(seqs)}] {time.time()-t0:.1f}s elapsed")
    print(f"  PromoGen2 embedding extraction: {len(seqs)} seqs, {n_params} params, "
          f"{time.time()-t0:.1f}s total, device={device}")
    return out


EXTRACTORS = {
    "dnabert2": extract_embeddings_dnabert2,
    "promogen2": extract_embeddings_promogen2,
}


def main(model_tag):
    extractor = EXTRACTORS[model_tag]

    oligo_ids, lib_seqs, fold = load_library_sequences()
    print(f"Library: {len(lib_seqs)} sequences")
    lib_emb = extractor(lib_seqs)
    np.savez_compressed(EMB_DIR / f"{model_tag}_library.npz",
                         oligo_ids=np.array(oligo_ids), embeddings=lib_emb, fold=fold)
    print(f"Wrote {EMB_DIR / f'{model_tag}_library.npz'} shape={lib_emb.shape}")

    rs241_ids, rs241_seqs = load_rs241_sequences()
    print(f"RS241: {len(rs241_seqs)} sequences")
    rs241_emb = extractor(rs241_seqs)
    np.savez_compressed(EMB_DIR / f"{model_tag}_rs241.npz",
                         ids=np.array(rs241_ids), embeddings=rs241_emb)
    print(f"Wrote {EMB_DIR / f'{model_tag}_rs241.npz'} shape={rs241_emb.shape}")


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "dnabert2"
    main(tag)
