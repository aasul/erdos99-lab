"""Build the static Pages site without running or changing the research optimizer."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '.site'
OUT.mkdir(exist_ok=True)
shutil.copytree(ROOT / 'desktop', OUT / 'desktop', dirs_exist_ok=True,
                ignore=shutil.ignore_patterns('*.desktop'))
shutil.copytree(ROOT / 'results/initial/figures', OUT / 'results/initial/figures', dirs_exist_ok=True)
(OUT / 'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="0;url=desktop/index.html"><title>Dot Lab</title><a href="desktop/index.html">Open Dot Lab</a></html>')
(OUT / '.nojekyll').touch()
print(f'Built site: {OUT}')
