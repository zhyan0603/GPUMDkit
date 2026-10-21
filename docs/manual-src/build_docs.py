"""Render the Markdown tutorials in docs/tutorials/{en,zh}/ as site-styled pages
(docs/htmls/{en,zh}/<file>.html) and write assets/js/gk-docs.js.

Called by build.py. Needs the `markdown` and `pymdownx` (pymdown-extensions) packages,
which MkDocs Material already installs. The Markdown files stay the single source.
"""
import html
import json
import re
from pathlib import Path

import markdown
from markdown.extensions.toc import slugify_unicode

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent
REPO_EDIT = "https://github.com/zhyan0603/GPUMDkit/edit/main/docs/tutorials"

GROUPS = {
    "start": ("Start", "开始"),
    "tools": ("Toolbox", "工具箱"),
    "guides": ("Guides", "指南"),
    "ref": ("Reference", "参考"),
    "dev": ("Contribute", "贡献"),
}
# slug, English file, Chinese file, group, parent slug (sub-pages are listed under their parent)
PAGES = [
    ("tutorials", "index", "index", "start", None),
    ("quick_start", "quick_start", "快速入门", "start", None),
    ("remote_web_console", "remote_web_console", "远程网页控制台", "start", None),
    ("simulation_and_postprocessing", "simulation_and_postprocessing", "模拟与后处理", "start", None),
    ("format_conversion", "format_conversion", "格式转换", "tools", None),
    ("structure_sampling", "structure_sampling", "结构采样", "tools", None),
    ("workflow_scripts", "workflow_scripts", "工作流脚本", "tools", None),
    ("calculator_scripts", "calculator_scripts", "计算器脚本", "tools", None),
    ("analyzer_scripts", "analyzer_scripts", "分析工具", "tools", None),
    ("plot_scripts", "plot_scripts", "绘图脚本", "tools", None),
    ("plot_nep", "plot_nep", "NEP训练绘图", "tools", "plot_scripts"),
    ("plot_md_structure", "plot_md_structure", "MD与结构绘图", "tools", "plot_scripts"),
    ("plot_diffusion", "plot_diffusion", "扩散绘图", "tools", "plot_scripts"),
    ("plot_heat_phonon", "plot_heat_phonon", "热输运与声子绘图", "tools", "plot_scripts"),
    ("polar_material_analysis", "polar_material_analysis", "极性材料分析", "tools", None),
    ("active_learning_workflow", "active_learning_workflow", "主动学习工作流", "guides", None),
    ("custom_commands", "custom_commands", "自定义命令", "guides", None),
    ("command_reference", "command_reference", "命令参考", "ref", None),
    ("contributing", "contributing_to_gpumdkit", "贡献指南", "dev", None),
]
STEMS = {"en": {e for _, e, *_ in PAGES}, "zh": {z for _, _, z, *_ in PAGES}}
STEM = {"en": 1, "zh": 2}  # index of the file stem inside a PAGES row

UI = {
    "en": {"manual": "Manual", "on_page": "On this page", "edit": "Edit this page", "prev": "Previous", "next": "Next", "lang": "en", "other": "zh"},
    "zh": {"manual": "手册", "on_page": "本页目录", "edit": "编辑此页", "prev": "上一页", "next": "下一页", "lang": "zh-CN", "other": "en"},
}


def make_md():
    return markdown.Markdown(extensions=[
        "admonition", "attr_list", "tables", "md_in_html", "pymdownx.details", "pymdownx.highlight", "pymdownx.superfences",
        "toc", "sane_lists",
    ], extension_configs={"pymdownx.highlight": {"use_pygments": False}, "toc": {"slugify": slugify_unicode, "toc_depth": "2-3"}})


def strip_emoji(text):
    return re.sub(r"^[^\w一-鿿]+", "", text).strip()


def read_page(lang, name):
    src = (DOCS / "tutorials" / lang / f"{name}.md").read_text(encoding="utf-8")
    m = re.match(r'\s*<div align="center">\s*<h1>(.*?)</h1>\s*<p[^>]*>(.*?)</p>\s*</div>', src, re.S)
    if m:
        title, desc = strip_emoji(re.sub(r"<[^>]+>", "", m.group(1))), m.group(2).strip()
        body = src[m.end():]
    else:  # plain "# Title" followed by a paragraph: that paragraph stays in the body and its first sentence is the lead
        m = re.match(r"\s*# ([^\n]+)\n+(.+?)(?:\n\n|$)", src, re.S)
        if not m:
            raise SystemExit(f"{lang}/{name}.md: expected a <div align=center> header or '# Title' and a paragraph")
        title, para = strip_emoji(m.group(1)), re.sub(r"\s+", " ", m.group(2))
        desc = re.split(r"(?<=[.。])\s", para, maxsplit=1)[0]
        body = src[m.start(2):]
    body = re.sub(r"^\s*←.*\n", "", body, count=1)  # "Back to …" lines are replaced by the breadcrumb
    # <!-- annotated-shot: name --> is replaced by the annotated screenshot block in manual-src/shots/name.html
    body = re.sub(r"<!-- annotated-shot: (\w+) -->", lambda m: "\n" + (HERE / "shots" / f"{m.group(1)}.html").read_text(encoding="utf-8") + "\n", body)
    md = make_md()
    out = md.convert(body)
    toc = []
    for t in md.toc_tokens:
        toc.append((t["id"], strip_emoji(html.unescape(re.sub(r"<[^>]+>", "", t["name"]))), 2))
        for c in t.get("children", []):
            toc.append((c["id"], strip_emoji(html.unescape(re.sub(r"<[^>]+>", "", c["name"]))), 3))
    return title, desc, out, toc


def polish(out, lang):
    # code blocks: caption with the language and a copy button
    btn = '<button class="copy-btn" data-copy-pre aria-label="Copy"><svg class="cp"><use href="#i-copy"/></svg><svg class="ok"><use href="#i-check"/></svg></button>'

    def code(m):
        cap = f'<div class="cap"><i></i><i></i><i></i><span>{m.group(1)}</span></div>' if m.group(1) else ""
        return f'<div class="code">{cap}<pre><code>{m.group(2)}</code></pre>{btn}</div>'
    out = re.sub(r'<pre(?: class="highlight")?><code(?: class="language-(\w+)")?>(.*?)</code></pre>', code, out, flags=re.S)
    out = re.sub(r"<table>", '<div class="tbl-wrap"><table class="tbl">', out)
    out = out.replace("</table>", "</table></div>")
    # links to other tutorials: foo.md#x -> foo.html#x ; index.md -> the manual hub
    def link(m):
        name, frag = m.group(1), m.group(2) or ""
        if name not in STEMS[lang]:
            raise SystemExit(f"unknown tutorial link: {name}.md ({lang})")
        return f'href="{name}.html{frag}"'
    out = re.sub(r'href="(?:\./)?([^"#/:]+)\.md(#[^"]*)?"', link, out)
    out = re.sub(r'<a href="(https?://[^"]+)"', r'<a href="\1" target="_blank" rel="noopener"', out)
    out = re.sub(r'<h([23]) id="([^"]+)">', r'<h\1 id="\2">', out)
    out = re.sub(r"<p>\s*(<a id=\"[^\"]*\"></a>\s*)+</p>", "", out)
    return out


def prefix_urls(chunk, p="../../"):
    return re.sub(r'(href|src)="(?!https?:|#|//|mailto:)([^"]+)"', lambda m: f'{m.group(1)}="{p}{m.group(2)}"', chunk)


def build_all(layout):
    head, rest = layout.split("<main>", 1)
    foot = rest.split("</main>", 1)[1]
    titles = {}
    pages = {}
    for lang, idx in (("en", 1), ("zh", 2)):
        for slug, e, z, group, parent in PAGES:
            titles[(lang, slug)], desc, body, toc = read_page(lang, (e, z)[idx - 1])
            pages[(lang, slug)] = (desc, polish(body, lang), toc)

    for lang in ("en", "zh"):
        ui = UI[lang]
        for n, (slug, e, z, group, parent) in enumerate(PAGES):
            desc, body, toc = pages[(lang, slug)]
            title = titles[(lang, slug)]
            fname_ = (e, z)[lang == "zh"]
            # sidebar
            side, last = '<aside class="man-side" aria-label="Pages">', None
            for row2 in PAGES:
                s2, g2, p2 = row2[0], row2[3], row2[4]
                if g2 != last:
                    side += f"<h6>{GROUPS[g2][0 if lang == 'en' else 1]}</h6>"
                    last = g2
                cls = ("on" if s2 == slug else "") + (" sub" if p2 else "")
                side += f'<a href="{row2[STEM[lang]]}.html" class="{cls.strip()}"><span>{html.escape(titles[(lang, s2)])}</span></a>'
            side += "</aside>"
            tocs = "".join(f'<a href="#{i}" class="l{lv}">{html.escape(t)}</a>' for i, t, lv in toc)
            rail = f'<nav class="doc-toc" aria-label="{ui["on_page"]}"><h6>{ui["on_page"]}</h6>{tocs}</nav>' if toc else "<nav class=\"doc-toc\"></nav>"
            prev = PAGES[n - 1] if n else None
            nxt = PAGES[n + 1] if n + 1 < len(PAGES) else None
            pn = '<div class="doc-pn">'
            pn += f'<a href="{prev[STEM[lang]]}.html"><small>{ui["prev"]}</small><b>{html.escape(titles[(lang, prev[0])])}</b></a>' if prev else "<span></span>"
            pn += f'<a class="next" href="{nxt[STEM[lang]]}.html"><small>{ui["next"]}</small><b>{html.escape(titles[(lang, nxt[0])])}</b></a>' if nxt else "<span></span>"
            pn += "</div>"
            crumb = f'<a href="../../manual.html">{ui["manual"]}</a><span>/</span>'
            if parent:
                prow = next(r for r in PAGES if r[0] == parent)
                crumb += f'<a href="{prow[STEM[lang]]}.html">{html.escape(titles[(lang, parent)])}</a><span>/</span>'
            else:
                crumb += f'<span>{GROUPS[group][0 if lang == "en" else 1]}</span><span>/</span>'
            fname = (e, z)[lang == "zh"]
            edit = f'<a class="doc-edit" href="{REPO_EDIT}/{lang}/{fname}.md" target="_blank" rel="noopener">{ui["edit"]}</a>'
            main = (f'<section class="man-hero doc-hero"><div class="container"><div class="crumb">{crumb}<span>{html.escape(title)}</span></div>'
                    f'<h1 class="h-section">{html.escape(title)}</h1><p class="lead">{desc}</p></div></section>\n'
                    f'<div class="container man doc">{side}<article class="man-body chap doc-body">{body}{pn}{edit}</article>{rail}</div>')
            page = head
            page = re.sub(r"<title>.*?</title>", f"<title>{html.escape(title)} — GPUMDkit</title>", page, count=1)
            plain = re.sub(r"<[^>]+>", "", desc)
            page = re.sub(r'<meta name="description" content="[^"]*">', f'<meta name="description" content="{html.escape(plain, quote=True)}">', page, count=1)
            page = re.sub(r'<link rel="canonical" href="[^"]*">', f'<link rel="canonical" href="https://gpumdkit.cn/htmls/{lang}/{fname_}.html">', page, count=1)
            page = page.replace('<body data-page="manual">', f'<body data-page="doc" data-root="../../" data-lang="{lang}" data-alt="../{ui["other"]}/{(e, z)[lang == "en"]}.html">')
            page = page.replace('<html lang="en" data-theme="light">', f'<html lang="{ui["lang"]}" data-theme="light">')
            # the page language wins over the stored preference
            page = page.replace('document.documentElement.classList.add("js");', f'document.documentElement.classList.add("js"); try {{ localStorage.setItem("gk-lang", "{lang}"); }} catch (e) {{}}', 1)
            page = prefix_urls(page)
            page = page.replace('href="../../manual.html" class="active"', 'href="../../manual.html" class="active"')
            tail = prefix_urls(foot).replace("gk-manual.js", "gk-doc.js")
            out = DOCS / "htmls" / lang
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{fname_}.html").write_text(page + "<main>\n" + main + "\n</main>" + tail, encoding="utf-8")

    # page list for the command palette, the "Detailed pages" rows and the site pages
    paths = {}
    items = ['  index: ["manual.html", "manual.html", "Manual", "使用手册"],']
    for slug, e, z, *_ in PAGES:
        paths[slug] = (f"htmls/en/{e}.html", f"htmls/zh/{z}.html")
        items.append(f'  {slug}: ["{paths[slug][0]}", "{paths[slug][1]}", {json.dumps(titles[("en", slug)], ensure_ascii=False)}, {json.dumps(titles[("zh", slug)], ensure_ascii=False)}],')
    (DOCS / "assets" / "js" / "gk-docs.js").write_text(
        "/* Generated by docs/manual-src/build.py from docs/tutorials. Do not edit by hand.\n   key -> [English page, Chinese page, English title, Chinese title], relative to the site root. */\n"
        "window.GK_DOCS = {\n" + "\n".join(items) + "\n};\n", encoding="utf-8")

    # docs/htmls/index.html kept for old links (e.g. the README): it opens the manual
    (DOCS / "htmls" / "index.html").write_text(
        '<!DOCTYPE html>\n<html><head><meta charset="utf-8"><title>GPUMDkit</title><meta http-equiv="refresh" content="0; url=../manual.html">'
        '<link rel="canonical" href="https://gpumdkit.cn/manual.html"></head><body><a href="../manual.html">GPUMDkit manual</a></body></html>\n', encoding="utf-8")

    # sitemap
    import datetime
    today = datetime.date.today().isoformat()
    urls = ["", "index_zh.html", "gallery.html", "publications.html", "manual.html"]
    urls += [f"htmls/en/{e}.html" for _, e, *_ in PAGES] + [f"htmls/zh/{z}.html" for _, _, z, *_ in PAGES]
    (DOCS / "sitemap.xml").write_text('<?xml version="1.0" encoding="UTF-8"?>\n<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n' + "".join(
        f"    <url>\n         <loc>{html.escape('https://gpumdkit.cn/' + u)}</loc>\n         <lastmod>{today}</lastmod>\n    </url>\n" for u in urls) + "</urlset>\n", encoding="utf-8")
    return len(PAGES) * 2, paths
