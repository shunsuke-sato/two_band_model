set term postscript eps enhanced color font "Arial, 20" dl 5

gap = 9.0/27.211
mu = 0.4
pvc = 0.5*sqrt(gap/mu)
eps_p(x) = x>gap ? 4.0*sqrt(2.0)*mu**1.5*pvc**2*sqrt(x-gap)/x**2 : 0
eps_np(x) = x>gap ? 4.0*(mu/gap)**1.5*pvc**2*sqrt(x**2-gap**2)/x : 0
ev = 27.211

set key left
set output "eps_i.eps"
p [5:15] "zeps.para.out" u ($1*ev):3 title "calc. (parabolic)" w l\
,eps_p(x/27.211) title "analytical. (parabolic)" w l\
,"zeps.out" u ($1*ev):3 title "calc. (non-parabolic)" w l\
,eps_np(x/27.211) title "analytical. (non-parabolic)" w l

unset output

