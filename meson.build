project('waves_cw','cpp')
deps=[]
#gsl_dep = dependency('gsl',fallback:['gsl','libgsl_dep'])
openmp_dep = dependency('openmp',language:'cpp')
#deps+=openmp_dep
executable('waves', 
'waves.cpp', 
'solver.cpp', 
'problem.cpp', 
dependencies:deps)
#tst = executable('test','test.cpp','solver.cpp','problem-simple.cpp')
#test('convergence test',tst)
