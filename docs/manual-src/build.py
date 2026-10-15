#!/usr/bin/env python3
"""Assemble docs/manual.html and the tutorial pages in docs/htmls/ (needs `markdown` + pymdown-extensions).

    python3 docs/manual-src/build.py

Edit the pieces, not manual.html:
  layout.html          page frame: <head>, navigation and footer (copied from gallery.html)
  hero.html            title block at the top of the page
  chapters/NN-*.html   one file per chapter; the number in the file name sets the order

Every chapter starts with <section class="chap" id="..." data-group="English|中文"> and
contains <span class="chap-no"></span> and one <h2>. The chapter number and the sidebar
are generated; a chapter inherits the sidebar group of the previous one if data-group is
the same, so only the group name needs to be repeated.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "manual.html"


def main():
    layout = (HERE / "layout.html").read_text(encoding="utf-8")
    from build_docs import build_all
    n_docs, paths = build_all(layout)
    chapters = []
    for n, f in enumerate(sorted((HERE / "chapters").glob("*.html")), 1):
        html = f.read_text(encoding="utf-8")
        m = re.match(r'\s*<section class="chap" id="([^"]+)" data-group="([^"|]+)\|([^"]+)"', html)
        h2 = re.search(r"<h2>(.*?)</h2>", html, re.S)
        if not (m and h2 and '<span class="chap-no"></span>' in html):
            raise SystemExit(f"{f.name}: needs <section class=\"chap\" id data-group>, <h2> and <span class=\"chap-no\"></span>")
        html = re.sub(r'href="[^"]*" data-doc="([a-z_]+)"', lambda m: f'href="{paths[m.group(1)][0]}" data-doc="{m.group(1)}"' if m.group(1) in paths else m.group(0), html)
        html = html.replace('<span class="chap-no"></span>', f'<span class="chap-no">{n:02d}</span>', 1)
        chapters.append((m.group(1), (m.group(2), m.group(3)), h2.group(1), n, html))

    side, last = '<aside class="man-side" aria-label="Manual contents">', None
    for cid, grp, title, n, _ in chapters:
        if grp != last:
            side += f'<h6><span class="l-en">{grp[0]}</span><span class="l-zh">{grp[1]}</span></h6>'
            last = grp
        side += f'<a href="#{cid}"><i>{n:02d}</i><span>{title}</span></a>'
    side += '</aside>'

    page = layout
    page = (page.replace("{{HERO}}", (HERE / "hero.html").read_text(encoding="utf-8").strip())
                .replace("{{SIDEBAR}}", side)
                .replace("{{CHAPTERS}}", "\n".join(c[4].strip() for c in chapters)))
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT} ({len(chapters)} chapters) and {n_docs} tutorial pages in docs/htmls/")


if __name__ == "__main__":
    main()
