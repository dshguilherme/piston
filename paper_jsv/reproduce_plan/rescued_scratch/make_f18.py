import csv, math
R=r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\results\r0_revision\fsi"
O=r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\paper_jsv\prototype\figs\data\f18_fsi.dat"
rows=list(csv.DictReader(open(R+r"\fsi_summary.csv")))
order=["c30_k1","c90_k1p5","c120_k1p5","c30_k3","c90_k7","c120_k7","c30_k18"]
V=["vac","dev","phys"]
out=["idx theta kR exvac exdev exphys sinvac sindev sinphys svac sdev sphys bvac bdev bphys"]
for i,c in enumerate(order,1):
    d={r["design"]:r for r in rows if r["cell"]==c}
    th=float(d["uniform"]["theta0_deg"]); kr=float(d["uniform"]["kR"])
    ex=[min(float(d["XU"]["dB_vs_uniform_"+v]),float(d["XS"]["dB_vs_uniform_"+v])) for v in V]
    sn=[float(d["SU"]["dB_vs_uniform_"+v]) for v in V]
    s=[float(d["SU"]["regret_"+v+"_dB"]) for v in V]
    b=[min(float(d["SU"]["regret_"+v+"_dB"]),float(d["SX"]["regret_"+v+"_dB"])) for v in V]
    out.append(" ".join([str(i),"%g"%th,"%.2f"%kr]+["%.2f"%x for x in ex+sn+s+b]))
    print(c, ["%.1f"%x for x in ex+sn+s+b])
open(O,"w").write("\n".join(out)+"\n")
# reopt
r=list(csv.DictReader(open(R+r"\reopt_eval_c30_k3.csv")))
for x in r: print(x["design"],x["variant"],"%.2f"%float(x["dB_vs_uniform_same_variant"]),"%.2f"%float(x["dB_vs_uniform_H8000"]))
P={(x["design"],x["variant"]):float(x["P_Fx"]) for x in r}
su=[x for x in rows if x["cell"]=="c30_k3" and x["design"]=="SU"][0]
print("SU vs reopt dev", 10*math.log10(float(su["P_dev"])/P[("reopt_dev_XU","dev")]))
print("SU vs reopt phys", 10*math.log10(float(su["P_phys"])/P[("reopt_phys_XU","phys")]))
# ranges
lo=[x for x in rows if float(x["kR"])<=3.01]
for key in ["KEa_over_KE_dev","KEa_over_KE_phys"]:
    v=[float(x[key]) for x in lo]; print(key,min(v),max(v))
v=[float(x["KEa_over_KE_phys"]) for x in rows if x["design"] in ("XU","XS","SU","SX")]; print("phys KE X/S all",min(v),max(v))
ex=[float(x["ratio_dev"]) for x in lo if x["design"] in ("XU","XS")]; print("ex dev",min(ex),max(ex))
ex=[float(x["ratio_phys"]) for x in lo if x["design"] in ("XU","XS")]; print("ex phys",sorted(ex))
o=[float(x["ratio_dev"]) for x in rows if x["design"] in ("uniform","SU","SX","EU","EX")]; print("other dev",min(o),max(o))
o=[float(x["ratio_phys"]) for x in rows if x["design"] in ("uniform","SU","SX","EU","EX")]; print("other phys",min(o),max(o))
a=[float(x["ratio_dev"]) for x in rows if x["cell"]=="c30_k18"]; print("anchor",min(a),max(a), [10*math.log10(t) for t in (min(a),max(a))])
k7=[(x["cell"],x["design"],float(x["dB_dev_over_vac"]),float(x["dB_phys_over_vac"])) for x in rows if x["cell"] in("c90_k7","c120_k7")]
print("k7 max dB dev",max(abs(t[2]) for t in k7),"phys",max(abs(t[3]) for t in k7))
bal=[abs(float(x["balance_resid_dev"])) for x in rows]; print("bal max",max(bal))
vf=[float(x["vel_frac_in_fluid_basis"]) for x in rows if x["design"] in ("XU","XS","SU","SX")]; print("velfrac",min(vf))
pr=[float(x["Prad_over_Pin_dev"]) for x in rows if x["design"] in ("XU","XS") and float(x["kR"])<=3.01]; print("PradPin exact lowkR",min(pr),max(pr))
sg=[(x["cell"],x["design"],float(x["sigma_vac"]),float(x["sigma_dev"]),float(x["sigma_phys"])) for x in lo if x["design"] in ("XU","XS","uniform")]
for t in sg: print(t)
