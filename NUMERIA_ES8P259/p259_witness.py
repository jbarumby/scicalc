"""ES8-P259 Gate A: exact-rational Fourier mass bound for the P258 odd witness.

Requires the authenticated P258_GL60.json from GitHub Actions run 44858303395.
This proves positivity of ONE vector after a banded-floor update, not SPD.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import argparse
import json

PINNED = "1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c"
L = Q("0.972955074527656652552676371721590")
T0 = 850
PERIODS = 23
BITS = 320

def require(ok, message):
    if not ok:
        raise ArithmeticError(message)

def pi_interval(steps=14):
    # Machin: pi=16 atan(1/5)-4 atan(1/239).
    # Alternating finite sums with even/odd end index bound atan.
    def bounds(n):
        def partial(last):
            return sum((Q((-1)**k, (2*k+1)*n**(2*k+1))
                        for k in range(last+1)), Q(0))
        return partial(2*steps+1), partial(2*steps)
    a,b=bounds(5)
    c,d=bounds(239)
    lo,hi=16*a-4*d,16*b-4*c
    require(3<lo<hi<Q(22,7), "bad pi enclosure")
    return lo,hi

def bessel_polys(last_degree):
    """j_n(z)=(P_n(1/z) sin(z)+Q_n(1/z) cos(z))/z."""
    p=[[1],[0,1]]
    q=[[0],[-1]]
    for n in range(1,last_degree):
        for a in (p,q):
            out=[0]*(len(a[n])+1)
            for j,v in enumerate(a[n]): out[j+1]+=(2*n+1)*v
            for j,v in enumerate(a[n-1]): out[j]-=v
            a.append(out)
    return p,q

def sqrt_lower(n, bits=BITS):
    return Q(isqrt(n << (2*bits)), 1 << bits)

def read_p258(path):
    raw=Path(path).read_bytes()
    require(sha256(raw).hexdigest()==PINNED, "P258 SHA-256 mismatch")
    j=json.loads(raw)
    require(j["schema"]=="NUMERIA.ES8P258.RAISED_AUXILIARY.1", "wrong schema")
    require(j["cutoff"]==504 and j["modes_per_sector"]==377, "wrong geometry")
    require(j["L_exact"]==str(L) and j["beta_raised_exact"]=="1/4", "wrong floor/width")
    require(j["G504_GL60_sha256"] is not None, "missing authenticated G504")
    w=j["sectors"]["odd"]["negative_witness"]
    coords=w["coordinates_exact_dyadic"]
    require(len(coords)==29 and coords[-1]==["1",0], "wrong witness")
    v=[Q(int(m))*Q(2)**int(e) / (2**39) for m,e in coords]
    oldlo,oldhi=map(Q,w["exact_rational_replay"]["rayleigh_interval_exact"])
    require(oldlo<oldhi<0, "old Rayleigh not strictly negative")
    return v,oldlo,oldhi

def poly_eval(coeff, y):
    value=Q(0)
    for z in reversed(coeff):
        value=value*y+z
    return value

def derivative_upper(coeff,err,y):
    return poly_eval([i*(abs(coeff[i])+err[i])
                      for i in range(1,len(coeff))], y)

def amplitudes(v, bits=BITS):
    p,q=bessel_polys(2*len(v)-1)
    length=2*len(v)
    a,b=[Q(0)]*length,[Q(0)]*length
    ea,eb=[Q(0)]*length,[Q(0)]*length
    root_error=Q(1,2**bits)
    for i,x in enumerate(v):
        n=2*i+1
        normroot=sqrt_lower(2*n+1,bits)
        phase=x*(-1)**i
        for j,z in enumerate(p[n]):
            a[j]+=phase*normroot*z
            ea[j]+=abs(phase*z)*root_error
        for j,z in enumerate(q[n]):
            b[j]+=phase*normroot*z
            eb[j]+=abs(phase*z)*root_error
    return a,b,ea,eb

def floor_gate():
    # P250 and P258 certify beta>0, and Binet:
    # Psi4(t)>=beta+log(t/504)+1/504-1/t for t>=504.
    # log(850/504)>2*(850-504)/(850+504)=346/677>1/2.
    lower=Q(2*(850-504),850+504)
    require(lower>Q(1,2),"second floor not established")
    return str(lower)

def prove_highpass(v, count=PERIODS):
    a,b,ea,eb=amplitudes(v)
    plo,phi=pi_interval()
    pmid=(plo+phi)/2
    norm2=sum((z*z for z in v),Q(0))
    require(norm2>Q(1,100), "witness norm too small")
    lower=Q(0)
    details=[]
    for k in range(count):
        # The kth full period has endpoints 850+2*pi*k/L etc.
        lo_t=Q(T0)+2*k*plo/L
        hi_t=Q(T0)+2*(k+1)*phi/L
        ylo=1/(L*hi_t)
        yhi=1/(L*lo_t)
        y0=1/(L*T0+(2*k+1)*pmid)
        require(ylo<y0<yhi,"invalid period interval")
        dy=max(yhi-y0,y0-ylo)
        ak,bk=poly_eval(a,y0),poly_eval(b,y0)
        epsa=poly_eval(ea,y0)+dy*derivative_upper(a,ea,yhi)
        epsb=poly_eval(eb,y0)+dy*derivative_upper(b,eb,yhi)
        rad2=ak*ak+bk*bk
        require(1024*(epsa+epsb)**2<rad2,
                "amplitudes too uncertain in period "+str(k))
        # L2 norm of sin/cos on a full 2*pi/L period is sqrt(pi/L)*sqrt(rad2).
        # Error L2 is <sqrt(2*pi/L)*sqrt(rad2)/32.
        # Since sqrt(2)<2, squared energy >=(15/16)^2 * (pi/L)*rad2.
        # Plancherel and t<=hi_t give the rational lower:
        term=Q(225,128)*rad2/(L*L*hi_t*hi_t*norm2)
        lower+=term
        details.append({"index":k,
                        "lower_mass_contribution_decimal":f"{float(term):.11e}",
                        "amplitude_variation_ratio_squared_upper_decimal":
                            f"{float((epsa+epsb)**2/rad2):.11e}"})
    require(hi_t<1000,"period integration extends beyond 1000")
    return lower,details,hi_t

def audit(path):
    floor=floor_gate()
    v,oldlo,oldhi=read_p258(path)
    mass,details,last=prove_highpass(v)
    raylo=oldlo+mass/4
    require(raylo>0,"banded-floor update inconclusive on the P258 witness")
    # Exact rational positivity is established above; the values below are
    # conservative decimal-grid enclosures, not float-based sign tests.
    grid=lambda x: Q((x*10**28).numerator//(x*10**28).denominator,10**28)
    return {
      "schema":"NUMERIA.ES8P259.ODD_WITNESS_FLOOR_GATE_A.1",
      "input_sha256":PINNED,
      "floor_after_850":"1/2",
      "floor_logarithm_lower_exact":floor,
      "old_rayleigh_interval_display":[f"{float(oldlo):.17e}",f"{float(oldhi):.17e}"],
      "highpass_strict_lower_grid_exact":str(grid(mass)),
      "updated_rayleigh_strict_lower_grid_exact":str(grid(raylo)),
      "highpass_lower_display":f"{float(mass):.15e}",
      "updated_rayleigh_lower_display":f"{float(raylo):.15e}",
      "period_count":len(details),
      "last_frequency_upper_display":f"{float(last):.10f}",
      "witness_nonzero_modes":29,
      "single_witness_strictly_positive_after_update":True,
      "full_377_by_377_sector_certified_positive":False,
      "Weil_crossing_established":False,
      "Riemann_hypothesis_proved":False,
      "proof_obligations":"P258 exact old Rayleigh, Binet/P250 beta>0, "
                          "spherical Bessel recurrence, Machin atan bounds, Plancherel",
      "periods":details,
    }

def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--p258",type=Path,default=Path("replay/P258_GL60.json"))
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    result=audit(args.p258)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("P259 Gate A",result["updated_rayleigh_lower_display"])

if __name__=="__main__":
    main()
