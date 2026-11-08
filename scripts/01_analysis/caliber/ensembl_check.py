import os
# -*- coding: utf-8 -*-
import urllib.request, json

def check(doi):
    url = "https://api.crossref.org/works/" + doi
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "refcheck/1.0 (mailto:r@example.com)"})
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.load(r)["message"]
        title = (d.get("title") or [""])[0]
        journal = (d.get("container-title") or [""])[0]
        year = None
        for k in ("published-print","published-online","issued"):
            if k in d and d[k].get("date-parts"):
                year = d[k]["date-parts"][0][0]; break
        vol = d.get("volume",""); issue = d.get("issue",""); page = d.get("page","")
        authors = d.get("author",[])
        vauth = []
        for a in authors:
            fam = a.get("family",""); given = a.get("given","")
            inits = "".join(p[0] for p in given.replace("-"," ").split() if p)
            vauth.append(f"{fam} {inits}")
        astr = ", ".join(vauth[:6])
        if len(vauth)>6: astr += ", et al."
        return f"OK | {astr} | {title[:70]} | {journal} {year};{vol}({issue}):{page}"
    except Exception as e:
        return f"ERR {type(e).__name__}: {e}"

out = []
# Ensembl 2022 数据库引用
out.append("[ensembl2022] " + check("10.1093/nar/gkab1049"))
# 备选：Ensembl 2020
out.append("[ensembl2020] " + check("10.1093/nar/gkz966"))

open(os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "ensembl_check.txt"),"w",encoding="utf-8").write("\n".join(out))
print("done")
