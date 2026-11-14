import sys, numpy as np, pandas as pd
src = open("code/replicate_doyle.py").read().split("# ---- Table 1")[0]
src = src.replace('Path(__file__).resolve().parent.parent', 'Path.cwd()')
exec(src)
PUBA = {(6,"ldl"):5.86,(6,"hdl"):-0.48,(6,"total_chol"):4.94,(6,"trigly"):3.63,(6,"weight"):1.95,(6,"sys"):-0.38,(6,"dias"):-0.25,
        (12,"ldl"):3.75,(12,"hdl"):-0.33,(12,"total_chol"):-0.51,(12,"trigly"):-25.63,(12,"weight"):1.88,(12,"sys"):-2.49,(12,"dias"):0.81}
rows=[]
for (month,var),pub in PUBA.items():
    col=f"{var}_val_{month}mo"; dd=raw.dropna(subset=[col])
    p=long_panel(dd,var,month)
    xs=["treat_now","post_period","post_treat_now","post_baseline","time",*PRESPEC]
    res=[]
    for s in range(40):
        est,ses=[],[]
        for q in mvn_impute(p,["panel","time"],seed=s):
            q["post_baseline"]=q.post_period*q.groupby("subject_id").panel.transform("first")
            r,df=xtreg_fe(q,"panel",xs); est.append(r["post_treat_now"][0]); ses.append(r["post_treat_now"][1])
        res.append(rubin(est,ses,df)[:2])
    e=np.array(res)
    rows.append(dict(month=month,var=var,pub=pub,mi_lo=np.percentile(e[:,0],2.5),mi_med=np.median(e[:,0]),mi_hi=np.percentile(e[:,0],97.5),
                     se_med=np.median(e[:,1]),z_from_med=(pub-np.median(e[:,0]))/np.median(e[:,1]),inside=np.percentile(e[:,0],2.5)<=pub<=np.percentile(e[:,0],97.5)))
out=pd.DataFrame(rows); print(out.round(2).to_string()); out.to_csv(sys.argv[1],index=False)
