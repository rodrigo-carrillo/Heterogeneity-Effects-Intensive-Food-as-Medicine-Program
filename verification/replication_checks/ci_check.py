import pandas as pd, numpy as np, re, subprocess, sys
txt = subprocess.run(["pdftotext","-layout","literature/trial/jamainternal_doyle_2023_oi_230084_1706307118.54674.pdf","-"],capture_output=True,text=True).stdout.replace("−","-")
L = txt.splitlines()
num = r"(-?[\d ]*\d(?:\.\d+)?)"
# Table 3
t3 = pd.read_csv("outputs/replication/table3.csv")
names = ["hba1c","ldl","hdl","total_chol","trigly","weight","sys","dias"]
labs = ["HbA1c level","Low-density","High-density","Total (n","Triglycerides","Weight, kg","Systolic","Diastolic"]
start = [i for i,l in enumerate(L) if "Table 3. Effects of Program" in l][0]
blk = L[start:start+30]
rows=[]
m6 = [i for i,l in enumerate(blk) if l.strip()=="6 mo"][0]; m12=[i for i,l in enumerate(blk) if l.strip()=="12 mo"][0]
for month, seg in ((6, blk[m6:m12]), (12, blk[m12:])):
    for v,lab in zip(names,labs):
        line=[l for l in seg if lab in l][0]
        vals=[float(x) for x in re.findall(r"-?\d+\.\d+", line.split(")",1)[1])]
        r=t3[(t3.month==month)&(t3.outcome.str.startswith(v))]
        rk = r[r.outcome=="weight (lb)"].iloc[0] if v=="weight" else r.iloc[0]
        rkg= r[r.outcome=="weight (kg)"].iloc[0] if v=="weight" else r.iloc[0]
        mine=[rk.ctl_mean, rk.ctl_sd, rk.trt_mean, rk.trt_sd, rk.unadj_est, rk.unadj_lo, rk.unadj_hi, rkg.adj_est, rkg.adj_lo, rkg.adj_hi]
        if v=="weight":   # published means in true kg, SDs in the 0.45-kg column
            mine[0]/=2.20462; mine[2]/=2.20462; mine[1]=rkg.ctl_sd; mine[3]=rkg.trt_sd
        for k,(pv,mv) in enumerate(zip(vals,mine)):
            rows.append(dict(table="T3",month=month,outcome=v,stat=["ctl_mean","ctl_sd","trt_mean","trt_sd","unadj","unadj_lo","unadj_hi","adj","adj_lo","adj_hi"][k],pub=pv,rep=mv))
# Table 4 (hand-transcribed)
t4 = pd.read_csv("outputs/replication/table4.csv")
T4PUB = """
6 ed_or_ip_i -7.13 -15.30 1.04 | 6 op_i 5.23 0 10.45 | 6 ed_or_ip_n -0.04 -0.27 0.20 | 6 op_n 1.00 0.28 1.72
6 ed_or_ip_payment -4 -1614 1607 | 6 op_payment 1376 -45 2796 | 6 nf_payment -86 -341 170 | 6 oth_payment -41 -552 470
6 total_payment 1048 -1568 3664 | 6 drug_insulin 3.01 -3.69 9.70 | 6 drug_metformin 7.18 0.71 13.64 | 6 drug_sulf -0.75 -4.84 3.33
6 drug_glp1 12.77 6.47 20.07 | 6 drug_other 3.24 -3.17 9.66 | 6 drug_any_db 2.73 -0.89 6.36
12 ed_or_ip_i -7.56 -16.84 0.72 | 12 op_i 0.47 -3.14 4.08 | 12 ed_or_ip_n -0.16 -0.59 0.27 | 12 op_n 0.20 -1.01 1.41
12 ed_or_ip_payment -1332 -3957 1294 | 12 op_payment 1348 -486 3181 | 12 nf_payment -150 -607 307 | 12 oth_payment -211 -1015 592
12 total_payment -424 -4465 3617 | 12 drug_insulin 2.99 -4.47 10.44 | 12 drug_metformin 8.87 1.88 15.86 | 12 drug_sulf -0.84 -5.85 4.16
12 drug_glp1 8.73 0.86 1.66 | 12 drug_other 0.44 -6.45 7.32 | 12 drug_any_db 3.29 -0.64 7.21"""
for item in re.split(r"[|\n]", T4PUB):
    if not item.strip(): continue
    mo,v,e,lo,hi = item.split(); mo=int(mo)
    r=t4[(t4.month==mo)&(t4.outcome==v)].iloc[0]
    for s,pv,mv in (("adj",float(e),r.est),("adj_lo",float(lo),r.lo),("adj_hi",float(hi),r.hi)):
        rows.append(dict(table="T4",month=mo,outcome=v,stat=s,pub=pv,rep=mv))
d=pd.DataFrame(rows)
def decs(x):
    s=f"{x:.6f}".rstrip("0"); return len(s.split(".")[1]) if "." in s and s.split(".")[1] else 0
d["dec"]=d.pub.map(decs)
d["rounded"]=[round(r,dd) for r,dd in zip(d.rep,d.dec)]
d["match"]=np.where(np.isclose(d.rounded,d.pub,atol=1e-9),"exact",np.where((d.rep-d.pub).abs()<=10.0**-d.dec*1.01,"last digit","DIFF"))
mi=(d.table=="T3")&d.stat.str.startswith("adj")&(d.outcome!="hba1c")
d.loc[mi,"match"]="MI cell"
print(d.groupby(["table","match"]).size())
print(d[d.match.isin(["last digit","DIFF"])].to_string())
d.to_csv(sys.argv[1],index=False)
