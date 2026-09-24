import pymupdf, re, sys
d = pymupdf.open(sys.argv[1] if len(sys.argv)>1 else "main2b_pt.pdf")
def limpa(t):
    out=[]
    for l in t.split("\n"):
        l=l.strip()
        if not l or re.fullmatch(r"\d{1,4}", l): continue
        if "Preprint submitted" in l or l.startswith("Page "): continue
        if "Resposta da fotossíntese" in l or "Unidade amostral em" in l: continue
        out.append(l)
    return out
pg=[limpa(p.get_text()) for p in d]
ruim=0
for i in range(len(pg)-1):
    if not pg[i] or not pg[i+1]: continue
    fim=pg[i][-1]; ini=" ".join(pg[i+1][:2])
    if not re.search(r"[.:;!?]$", fim) and re.match(r"^(Table|Figure)\s*\d", ini):
        ruim+=1
        print(f"  p.{i+1}->{i+2}: '...{fim[-45:]}'  ||  '{ini[:40]}'")
print(f"paginas: {len(pg)} | cortes por float: {ruim}")
