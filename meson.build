project('waves_cw','cpp')
deps=[]
#gsl_dep = dependency('gsl',fallback:['gsl','libgsl_dep'])
openmp_dep = dependency('openmp',language:'cpp')
deps+=openmp_dep
executable('waves','waves.cpp','problem.cpp',dependencies:deps)
executable('simple','waves.cpp','problem-simple.cpp',dependencies:deps)
