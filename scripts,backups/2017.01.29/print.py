a = ["ERROR CODE", "start time (YYYYMMDD-hhmmssms)", "v0 (kpc/s)", "mag(v0) (kpc/s)", "mag(v0/c)", "E0 (kg*kpc^2/s^2)", "theta0 (rads)", "phi0 (rads)", "q_b (A*kpc)", "mass (kg)", "dt (s)", "distance to halt (kpc)", "end time (YYYYMMDD-hhmmssms)","vel_final (kpc/s)", "mag(vel_final) (kpc/s)", "mag(vel_final/c)", "pos_final(kpc)", "mag(pos_final) (kpc)", "theta_final (rads)", "phi_final (rads)", "distance tracked (kpc)", "distance from start (kpc)","KE_final (kg*kpc^2/s^2)", "max velocity (kpc/s)", "max velocity/c","arclength in rho1kpc region (kpc)", "time after t=0 start (s)","iterations", "real runtime (real s)", "acc_final (kpc/s/s)","exit status (True/False)"]

print()
for n in range(len(a)):
    if n<10:
        print("[%d]  %s" % (n,a[n]))
    else:
        print("[%d] %s" % (n,a[n]))