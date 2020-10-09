//Resolution complete:PUT
Complete_put=fscanfMat("C:\Users\imane\eclipse-workspace\bs\solution_complete_put.txt")
Complete_put(1)=[]
Complete_put(1)=[]
plot(linspace(0,300,length(Complete_put)),Complete_put,'LineWidth',3)
xgrid
xtitle("Approche de C(0,S):EDP complete ","x","Put" )
xs2png(0,'complete_put.png');
scf()
//Resolution complete:CALL
Complete_call=fscanfMat("C:\Users\imane\eclipse-workspace\bs\solution_complete_call.txt")
plot(linspace(0,300,length(Complete_call)),Complete_call,'LineWidth',3)
xgrid
xtitle("Approche de C(0,S):EDP complete " ,"x","Call" )
xs2png(0,'complete_call.png');
scf()
//Resolution reduite:PUT
Reduite_put=fscanfMat("C:\Users\imane\eclipse-workspace\bs\solution_reduite_put.txt")
plot(linspace(0,300,length(Reduite_put)),Reduite_put,'LineWidth',3)
xgrid
xtitle( "Approche de C(0,S):EDP réduite ","x","put" )
xs2png(0,'reduite_put.png');
scf()
//Resolution reduite:CALL
Reduite_call=fscanfMat("C:\Users\imane\eclipse-workspace\bs\solution_reduite_call.txt")
plot(linspace(0,300,length(Reduite_call)),Reduite_call,'LineWidth',3)
xgrid
xtitle( "Approche de C(0,S):EDP réduite ","x","call" )
xs2png(0,'reduite_call.png');


