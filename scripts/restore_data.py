#!/usr/bin/env python3
"""Restore and SHA-256 verify all competition, external, and historical scored files
using registry/data_manifest.json. Uses `gh api` for authentication.
"""
import argparse, hashlib, json, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def sha256_file(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def fetch(repo, ref, path, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + '.partial')
    try:
        with tmp.open('wb') as out:
            subprocess.run(
                ['gh','api', f'repos/{repo}/contents/{path}?ref={ref}',
                 '-H','Accept: application/vnd.github.raw'],
                stdout=out, check=True, timeout=300)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)

def restore_entry(entry, root):
    dest = root / entry['dest']
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and sha256_file(dest) == entry['sha256']:
        return {'id':entry['id'],'dest':entry['dest'],'group':entry['group'],
                'status':'present','bytes':dest.stat().st_size,'sha256':entry['sha256']}
    if 'parts' in entry:
        tmp = dest.with_suffix('.assembling')
        with tmp.open('wb') as out:
            for part in entry['parts']:
                p = root / 'raw_parts' / Path(part).name
                p.parent.mkdir(parents=True, exist_ok=True)
                fetch(entry['repo'], entry['ref'], part, p)
                with p.open('rb') as src:
                    shutil.copyfileobj(src, out, 1 << 20)
                p.unlink()
        got = sha256_file(tmp)
        if got != entry['sha256']:
            tmp.unlink(missing_ok=True)
            raise SystemExit(f'HASH MISMATCH for {entry["id"]}: got {got} expected {entry["sha256"]}')
        tmp.replace(dest)
    else:
        fetch(entry['repo'], entry['ref'], entry['path'], dest)
        got = sha256_file(dest)
        if got != entry['sha256']:
            dest.unlink(missing_ok=True)
            raise SystemExit(f'HASH MISMATCH for {entry["id"]}: got {got} expected {entry["sha256"]}')
    return {'id':entry['id'],'dest':entry['dest'],'group':entry['group'],
            'status':'restored','bytes':dest.stat().st_size,'sha256':entry['sha256']}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', default='core', choices=['core','external','scored','all'])
    ap.add_argument('--manifest', default=str(ROOT/'registry'/'data_manifest.json'))
    ap.add_argument('--target-dir', default=str(ROOT/'data'))
    args = ap.parse_args()
    man = json.loads(Path(args.manifest).read_text())
    root = Path(args.target_dir)
    root.mkdir(parents=True, exist_ok=True)
    results = []
    for e in man['files']:
        if args.group in ('all', e['group']):
            r = restore_entry(e, root)
            results.append(r)
            print(f"{r['status']:>8}  {e['dest']:60s}  {r['bytes']:>10d}B  sha={r['sha256'][:12]}")
    # Preserve verified results from earlier --group calls while checking that
    # each prior file still exists and matches its pinned hash. Write the receipt
    # beneath the selected storage root, not unconditionally under ROOT/data.
    receipt_path = root / 'restore_receipt.json'
    merged = {}
    if receipt_path.is_file():
        try:
            old = json.loads(receipt_path.read_text())
            for item in old.get('files', []):
                dest = root / item.get('dest', '')
                if dest.is_file() and sha256_file(dest) == item.get('sha256'):
                    merged[item.get('id', item.get('dest'))] = item
        except (OSError, ValueError, TypeError):
            merged = {}
    for item in results:
        merged[item['id']] = item
    all_files = sorted(merged.values(), key=lambda item: item.get('id', ''))
    receipt = {
        'schema_version': 2,
        'verified_files_count': len(all_files),
        'storage_root': str(root),
        'last_requested_group': args.group,
        'files': all_files,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'Verified {len(results)} requested files; {len(all_files)} verified files in {receipt_path}')

if __name__ == '__main__':
    main()
