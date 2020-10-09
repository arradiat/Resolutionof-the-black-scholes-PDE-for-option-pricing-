//Erreur:PUT
erreur_put=fscanfMat("C:\Users\imane\eclipse-workspace\bs\erreur_put.txt")
plot(linspace(0,300,length(erreur_put)),erreur_put,'LineWidth',3)
xgrid
xtitle( "Erreur:Put ","x","Erreur" )
xs2png(0,'erreur_put.png');
scf()
//Erreur:CALL
erreur_call=fscanfMat("C:\Users\imane\eclipse-workspace\bs\erreur_call.txt")
plot(linspace(0,300,length(erreur_call)),erreur_call,'LineWidth',3)
xgrid
xtitle( "Erreur:Call","x","Erreur" )
xs2png(0,'erreur_call.png');
