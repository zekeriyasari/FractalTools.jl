using FractalTools
using Plots
ifs = Sierpinski()
initset = [rand(2)]

atr = attractor(ifs, initset, alg=DetAlg(), numiter=5) 
atr.ifs
atr.initset
atr.numiter
atr.set

scatter(getindex.(atr.set,1) , getindex.(atr.set,2))