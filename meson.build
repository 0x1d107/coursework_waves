project('waves_cw','cpp')
deps=[]
#gsl_dep = dependency('gsl',fallback:['gsl','libgsl_dep'])
openmp_dep = dependency('openmp',language:'cpp')
deps+=openmp_dep
executable('waves','waves.cpp','problem.cpp',dependencies:deps)
executable('simple','waves.cpp','problem-simple.cpp',dependencies:deps)
executable('analytic','waves.cpp','problem-analytic.cpp',dependencies:deps)
executable('refraction','waves.cpp','problem-refraction.cpp',dependencies:deps)
executable('difraction','waves.cpp','problem-difraction.cpp',dependencies:deps)
executable('detector','waves.cpp','problem-detector.cpp',dependencies:deps)
