import sys, numpy as np, pandas as pd
src = open("code/replicate_doyle.py").read().split("# ---- Table 1")[0].replace('Path(__file__).resolve().parent.parent','Path.cwd()')
exec(src)
PUBA = {(6,"ldl"):5.86,(6,"hdl"):-0.48,(6,"total_chol"):4.94,(6,"trigly"):3.63,(6,"weight"):1.95,(6,"sys"):-0.38,(6,"dias"):-0.25,
        (12,"ldl"):3.75,(12,"hdl"):-0.33,(12,"total_chol"):-0.51,(12,"trigly"):-25.63,(12,"weight"):1.88,(12,"sys"):-2.49,(12,"dias"):0.81}
xs=["treat_now","post_period","post_treat_now","post_baseline","time",*PRESPEC]
rows=[]
for (month,var),pub in PUBA.items():
    dd=raw.dropna(subset=[f"{var}_val_{month}mo"]); p=long_panel(dd,var,month)
    miss_ids=p.loc[(p.post_period==0)&p.panel.isna(),"subject_id"].values
    for variant,kmin in (("A: one draw reused",1),("B: draw or missing",0)):
        res=[]
        for s in range(40):
            rng=np.random.default_rng(1000+s); imps=mvn_impute(p,["panel","time"],seed=s)
            base0=p[p.post_period==0].set_index("subject_id").panel
            pick={i:rng.integers(kmin,11) for i in miss_ids}
            chosen=base0.copy()
            for i,k in pick.items():
                chosen[i]=np.nan if k==0 else imps[k-1].query("post_period==0").set_index("subject_id").panel[i]
            est,ses=[],[]
            for q in imps:
                q=q.copy(); q["post_baseline"]=np.where(q.post_period.eq(1),q.subject_id.map(chosen),0.0)
                q=q[~(q.post_period.eq(1)&q.post_baseline.isna())]
                r,df=xtreg_fe(q,"panel",xs); est.append(r["post_treat_now"][0]); ses.append(r["post_treat_now"][1])
            res.append(rubin(est,ses,df)[0])
        e=np.array(res)
        rows.append(dict(month=month,var=var,variant=variant,pub=pub,lo=np.percentile(e,2.5),med=np.median(e),hi=np.percentile(e,97.5),
                         inside=np.percentile(e,2.5)<=pub<=np.percentile(e,97.5)))
out=pd.DataFrame(rows); print(out.round(2).to_string()); out.to_csv(sys.argv[1],index=False)
