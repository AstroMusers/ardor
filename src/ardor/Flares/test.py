import ardor.Flares.Flare as Flare
import ardor.Flares.allesfitter_priors as priors

a, b,c = priors.flare_energy(0.001, 0.1, 5550,0.8, uncertainty = True)

print(a,b,c)
