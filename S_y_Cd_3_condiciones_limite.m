%FILTROS: TIPOS DE PRUEBAS
clear;clc;

Tiempo_produccion_en_PBU = 72; % Solo se usa este valor si FILTRO2 = 2
FILTRO1 = 1; %  1: petroleo, 2: GAS
FILTRO2 = 2; %  1: Prueba DD, 2: prueba PBU
FILTRO3 = 1; %  2: Comportamiento infinito, 2: Presion constante, 3: Barrera de no flujo,...
             %  4: Com. inf. yacimiento con fracturas, 5: Fluido no
             %  newtoniano inf, 6: Fluido no newtoniano presion constante,
             %  7: fluido no newtoniano barrera de no flujo
FILTRO4 = 1; %  1: Yacimiento no sensitivo, 2: yacimiento sensitivo, 3: Sensitivo y no sensitivo
FILTRO5 = 2; %  1: Ajuste con Cd, 2: sin ajuste de Cd
FILTRO6 = 1; %  1: Ajuste con K, 2 : Sin ajuste de K
FILTRO7 = 2; %  1: Ajuste con S, 2 : Sin ajuste de S
FILTRO8 = 4; %  1: gráfico cartesiano pw vs t, 2 : Curva de la derivada
             %  3 : grafico semilog, 4 : Curva de la derivada & Rammey
FILTRO9 = 3; %  1: Resultados reales, 2: Resultados analíticos, 3: Resultados reales y analíticos 
%%
%FILTRO: NUMERO DE LA PRUEBA (DATOS)
PBU=3; % Número de la prueba (1:4)PBU Oil, (4:5) PBU_gas; escoja 0 si no se utilizará
DD=0; %Número de la prueba (1:3) escoja 0 si no se utilizará
Matriz=DATOS_PT(PBU,DD);
datosyac=DATOS_PP(PBU,DD);
%%
%CREACIÓN DE MATRICES
vector_k = (1:1:10);%22 mD
vector_S = (0:2:20);
vector_Cd =(1:1:3);
vector_Cd = 10.^vector_Cd;
vector_gamma = (-6:0.5:-3);
vector_gamma = 10.^vector_gamma;
vector_w =(-2:0.5:-4);
vector_w = 10.^vector_w;
vector_lambda = (-4:0.5:-6);
vector_lambda = 10.^vector_lambda;
%% DATOS_ROCA
P_i = datosyac(1,1);
phi = datosyac(2,1);
C = datosyac(3,1);
API=datosyac(4,1);
r_w = datosyac(5,1);
q = datosyac(6,1); 
h = datosyac(7,1);
r_e=datosyac(8,1);
tp=datosyac(9,1);
Temp=datosyac(10,1); 
Red=r_e/r_w;
Temp=Temp + 459.67;
nu=0;%Parámetro del fluido no newtoniano
%%
%Calculo de beta y v del fluido no newtoniano
v=(1-nu)/2;
beta=(2*v)/(3-nu);
%Calculo Red
if FILTRO1 == 1 %OIL
 Pb=900;
end
gammag_oil=0.65;%GAS ASOCIADO AL PETRoLEO
Gamma_C=141.5/(131.5+API);
vector_t = Matriz(1:length(Matriz),1);
vector_presion = Matriz(1:length(Matriz),2);
longitud_vector_k = length(vector_k);
longitud_vector_t = length(vector_t);
longitud_vector_S = length(vector_S);
longitud_vector_Cd = length(vector_Cd);
longitud_vector_gamma = length(vector_gamma);
longitud_vector_w = length(vector_w);
longitud_vector_lambda= length(vector_lambda);
vector_presion_calculada = zeros(longitud_vector_t,1);
P_inicial = P_i;
deltap = zeros(longitud_vector_t,1);
if FILTRO1 == 2 %GAS
    vector_MP = zeros(longitud_vector_t,1);
    vector_presion_interpolacion = (14.7:50:2*P_i);
    for i=1:length(vector_presion_interpolacion)
        vector_MP(i) = MP(Gamma_C,vector_presion_interpolacion(i), Temp,500);
    end
    P_i = MP(Gamma_C,P_i,Temp,500);
end
%%
%AJUSTE DE CD
Ct=6.16E-05;
Bgas = 0.7364;
dtedp=60;
tDdpD=7;
Qref= 6939.68;
Bbifasico = 0.7364;
mugas=0.0771;
Pasterisco=5700.1;
Kgas=8  ;
if FILTRO5 == 1
    Presion1 = Matriz(4,2);
    Presion2 = Matriz(3,2);
    Tiempo1  = Matriz(4,1);
    if FILTRO2 == 2
        deltate1 = Tiempo1*tp/(Tiempo1+tp);
        deltape1 = Presion1-Presion2;
    end
     if FILTRO2 == 1
        deltate1 = Tiempo1;
        deltape1 = Presion1-Presion2;
    end
    CD = (0.037233*Qref*Bgas)/(h*phi*Ct*r_w^2)*deltate1/deltape1;    
end
%AJUSTE VALOR DE K
%Punto de match
if FILTRO6 == 1
    K = 162.6*Qref*mugas*Bbifasico/h*(tDdpD/dtedp);    
end
%Ajuste del daño S
if FILTRO7 == 1
    Ss=Kgas*h*(Pasterisco-P_inicial)/(141.2*Qref*mugas*Bgas);
end
%%
tic
if FILTRO4 == 1 %Yacimiento no sensitivo
    l_iteracion = 1;
    m_iteracion = 1;
elseif FILTRO4 == 2 % Yacimiento sensitivo
    l_iteracion = 2;
    m_iteracion = 2;
else % No sensitivo y sensitivo
    l_iteracion = 2;
    m_iteracion = 1;
end

if FILTRO4== 2 || FILTRO4== 3 
    KIN = 186; % Solo en prueba sensitiva
    S_total = 6; % Solo en prueba sensitiva
end
for sen=m_iteracion:l_iteracion
    tic
    if sen == 1
        SENSITIVO = 0;
        longitud_k = longitud_vector_k;
    else
        SENSITIVO = 1;
        longitud_k = longitud_vector_gamma;
    end
    if SENSITIVO == 1 %#ok<*STCMP>
        if FILTRO4 == 3
        KIN  = PERMEABILIDAD;
        end
        if FILTRO4 == 3
        S_lim=S_total; %S_total= SM + SG, si no sensitivo SG=0;
        vector_S = (0:1:S_lim);%S ANTERIOR ES EL DAÑO TOTAL; DAÑO MÁXIMO
        longitud_vector_S = length(vector_S);
        end
    end
    sigma = 1000;
    if FILTRO6==1
        longitud_k = 1;     
    end
    for k = 1:longitud_k
        if SENSITIVO == 0%no sensitivo
            if FILTRO6==1%con ajuste de K
               permeabilidad =  K;
            else
               permeabilidad = vector_k(k);
            end
        else
            gamma = vector_gamma(k);
        end
        if FILTRO7==1
          longitud_vector_S=1;     
        end
        for s = 1:longitud_vector_S
            if FILTRO7==1
                S=Ss;
            else
                S = vector_S(s);
            end
            
            if FILTRO5==1
                longitud_vector_Cd = 1;    
                Cd=CD; 
            end
            for c = 1:longitud_vector_Cd 
                if FILTRO5==2
                 Cd = vector_Cd(c);
                end
                for t = 1:longitud_vector_t
                    if FILTRO1 == 1 % OIL
                        mu = functionmuo(vector_presion(t),Temp,gammag_oil,API,Pb); % Revisar módulo PVT
                        B = functionBo(vector_presion(t),Temp,gammag_oil,API,Pb); % Revisar módulo PVT
                    end
                    if FILTRO1 == 2 % GAS
                        mu = functionmug(vector_presion(t),Temp,Gamma_C); % Revisar módulo PVT
                        B  = functionBg(vector_presion(t),Temp,Gamma_C); %BY/PCN % Solo en el caso de sensitivo
                    end
                    if SENSITIVO == 0
                         tD = 0.000264*permeabilidad*vector_t(t)/(phi*mu*C*r_w^2);
                    
                         if FILTRO3 == 1% Comportamiento infinito                     
                           Pd = steh('Solucion_comportamiento_infinito',tD,S,Cd,Red);  
                         end
                         if FILTRO3 == 4%inf. con fracturas 
                               for w = 1:longitud_vector_w
                                 W=vector_w(i);
                                 for l = 1:longitud_vector_lambda
                                   L=vector_lambda(l);
                                   Pd = steh_fracturas('Yac_fracturado_infinito',tD,S,Cd,Red,W,L);
                                 end
                               end  
                         end
                         if FILTRO3 == 5% inf. no newtoniano 
                             Pd = steh_no_newtoniano('sol_inf_no_newtoniano',tD,S,Cd,Red,v,beta);  
                         end
                         if FILTRO3 == 6% presion_constante_no newtoniano 
                             Pd = steh_no_newtoniano('sol_PC_no_newtoniano',tD,S,Cd,Red,v,beta);  
                         end
                         if FILTRO3 == 7% noflujo_no newtoniano 
                             Pd = steh_no_newtoniano('sol_noflujo_no_newtoniano',tD,S,Cd,Red,v,beta);  
                         end
                         if FILTRO3 == 2% Barrera de presion constante  
                           Pd = steh('Solucion_presion_constante',tD,S,Cd,Red);            
                         end
                         if FILTRO3 == 3% Barrera de no flujo    
                           Pd = steh('Solucion_barrera_no_flujo',tD,S,Cd,Red);  
                         end
                    
                           if FILTRO2 == 2 % PBU
                             tDtotal = 0.000264*permeabilidad*(vector_t(t)+tp)/(phi*mu*C*r_w^2);
                             if FILTRO3 ==1% Comportamiento infinito
                               Pd = steh('Solucion_comportamiento_infinito',tDtotal,S,Cd,Red) - Pd;
                             end
                             if FILTRO3 == 4%inf. con fracturas 
                                for w = 1:longitud_vector_w
                                    W=vector_w(i);
                                    for l = 1:longitud_vector_lambda
                                         L=vector_lambda(l);
                                         Pd = steh_fracturas('Yac_fracturado_infinito',tDtotal,S,Cd,Red,W,L) - Pd;
                                    end
                                end  
                             end
                             if FILTRO3 == 5%inf. no newtoniano 
                               Pd = steh_no_newtoniano('sol_inf_no_newtoniano',tDtotal,S,Cd,Red,v,beta) - Pd;  
                             end
                             if FILTRO3 == 6% presion_constante_no newtoniano 
                              Pd = steh_no_newtoniano('sol_PC_no_newtoniano',tDtotal,S,Cd,Red,v,beta)- Pd;  
                             end
                             if FILTRO3 == 7% noflujo_no newtoniano 
                              Pd = steh_no_newtoniano('sol_noflujo_no_newtoniano',tDtotal,S,Cd,Red,v,beta)- Pd;  
                             end
                             if FILTRO3 == 2% Barrera de presion constante
                                Pd = steh('Solucion_presion_constante',tDtotal,S,Cd,Red) - Pd;
                             end
                             if FILTRO3 == 3% Barrera de no flujo
                                 Pd = steh('Solucion_barrera_no_flujo',tDtotal,S,Cd,Red) - Pd;   
                             end
                           end 
                        
                           if FILTRO1 ==1 % OIL
                            vector_presion_calculada(t)= P_i - (141.2*Pd*q*mu*B/(permeabilidad*h));%P(t)
                           end
                           if FILTRO1 == 2 % GAS, Pseudopresion calculada: Pi convertido a valor inicial de mP(i)
                            vector_presion_calculada(t)= P_i - (1422*q*Temp*Pd /(permeabilidad*h) );%mP(t)
                           end
                    end
                    if SENSITIVO == 1
                        permeabilidad = KIN * exp( -gamma*(P_inicial- vector_presion(t)) );
                        tD = 0.000264*permeabilidad*vector_t(t)/(phi*mu*C*(r_w^2));

                        if FILTRO1 == 1
                            gammaD = 141.2*gamma *q*mu*B / (KIN*h);
                        end
                        if FILTRO1 == 2
                            gammaD = 141.2*gamma *q*mu*B*1000 /(KIN*h);
                        end
                        if FILTRO2 == 1 %DDN
                            if FILTRO3 == 1% Comportamiento infinito                     
                             UO1S = steh('Solucion_comportamiento_infinito',tD,S,Cd,Red);   
                             UOWPROM=UO1S;
                            end
                            if FILTRO3 == 2% Barrera de presion constante  
                             UO1S = steh('Solucion_presion_constante',tD,S,Cd,Red);  
                             UOWPROM=UO1S;
                            end
                            if FILTRO3 == 3% Barrera de no flujo    
                             UO1S = steh('Solucion_barrera_no_flujo',tD,S,Cd,Red); 
                             UOWPROM=UO1S;
                            end
                        else %PBU

                            Ttotal = tp + vector_t(t);
                            DTD= (2.64E-4)* Ttotal * permeabilidad /( phi*mu*C*(r_w^2) );%Tp+delta t adimen.
                            if FILTRO3 ==1% Comportamiento infinito
                               UO1S2 = steh('Solucion_comportamiento_infinito',tD,S,Cd,Red);
                               UO1S1 = steh('Solucion_comportamiento_infinito',DTD,S,Cd,Red);
                               UOWPROM=  UO1S1 - UO1S2;
                             end
                             if FILTRO3 == 2% Barrera de presion constante
                               UO1S2 = steh('Solucion_presion_constante',tD,S,Cd,Red);
                               UO1S1 = steh('Solucion_presion_constante',DTD,S,Cd,Red);
                               UOWPROM=  UO1S1 - UO1S2;
                             end
                             if FILTRO3 == 3% Barrera de no flujo
                               UO1S2 = steh('Solucion_barrera_no_flujo',tD,S,Cd,Red);
                               UO1S1 = steh('Solucion_barrera_no_flujo',DTD,S,Cd,Red);  
                               UOWPROM=  UO1S1 - UO1S2;
                             end   
                        end
                        if (gammaD*UOWPROM) < 1
                            PWFD = (-1/gammaD)*log( 1 - ( gammaD*UOWPROM ) ); % Esta si es PD
                        else
                            PWFD=NaN;
                        end
                        if FILTRO1 == 1
                            vector_presion_calculada(t)= P_i - (141.2*q*mu*B*PWFD /(permeabilidad*h) );%P(t)
                        end
                        if FILTRO1 == 2
                            vector_presion_calculada(t)= P_i - (1422*q*Temp*PWFD /(permeabilidad*h) );%mP(t)
                        end
                    end
                    if FILTRO1 == 2
                        vector_presion_calculada(t) = interp1(vector_MP,vector_presion_interpolacion,vector_presion_calculada(t));
                    end
                end
                sigma_calculado = sqrt( (sum( ((vector_presion - vector_presion_calculada)./vector_presion).^2 ))...
                    /longitud_vector_t )*100
                xx = [k ,s , c]%iteraci�n de cada variable
                if(sigma_calculado < sigma)
                    if SENSITIVO == 0
                        PERMEABILIDAD = permeabilidad;
                        S_total = S;
                        sigma_sin_sensitivo = sigma_calculado;
                        CsD_sin_sensitivo = Cd;
                        CsD=CsD_sin_sensitivo;
                        Presion_sin_sensitivo = vector_presion_calculada;
                        Presion_analitica = Presion_sin_sensitivo;
                        S_lim = S_total; % Este valor solo servira si FILTRO3 = 3
                    else
                        GAMMA = gamma;
                        SM = S; %Daño Mecánico
                        sigma_con_sensitivo = sigma_calculado;
                        Presion_con_sensitivo = vector_presion_calculada;
                        Presion_analitica = Presion_con_sensitivo;
                        CsD_con_sensitivo = Cd;
                    end
                    sigma = sigma_calculado;
                end
            end
        end
    end
    %%
    if SENSITIVO == 1%sensitivo
        SG=S_total - SM;
        if FILTRO1 == 2
            P_i = P_inicial;
        end
        K_con_sensitivo = KIN*exp( -GAMMA*(P_i - vector_presion(longitud_vector_t) ) );
    end
    toc
end

%%
disp('TIEMPO TOTAL DE ITERACIoN');
toc
disp('  ');
if FILTRO4 == 1
 K = PERMEABILIDAD;
 fprintf('Permeabilidad=%2.0f mD \n',K);
 fprintf('CsD=%2.0f \n',CsD);
 fprintf('S_total=%2.1f \n',S_total);
 fprintf('Sigma =%2.2f \n \n',sigma);
 fprintf('bsd=%2.0f  \n',xx);
end
vector_Presion = Matriz(:,2);
vector_t = Matriz(:,1);
data=[vector_t' vector_Presion'];
if FILTRO4 == 2 || FILTRO4 == 3
 disp('SENSITIVO');
 fprintf('Gamma= %4.1e \n',GAMMA);fprintf('K =%2.0f mD \n',K_con_sensitivo);
 fprintf('CsD=%2.0f \n',CsD_con_sensitivo);
 fprintf('S_Mecánico=%2.1f \t',SM);fprintf('S_Geomecánico=%2.1f \n',SG);
 fprintf('Sigma =%2.2f \n',sigma_con_sensitivo);
end
%%
%Vector t para graficar la derivada
vector_T=zeros(longitud_vector_t-1,1);
for i=1:longitud_vector_t-1
    vector_T(i,1) = vector_t(i+1,1);
end
%Gráfico pw vs t
if FILTRO8 ==1
    plot1=plot(vector_t,vector_Presion,'bo','LineWidth', 1.5);
    hold on
    plot2 = plot(vector_t,Presion_analitica,'rs','LineWidth',1.5);
    hold on
    title( 'Prueba de Presion' , 'FontSize' ,12 )
    grid on
    xlabel ( 't [h] ')
    ylabel ( 'Pw [lpc]')
    legend('P.real','P.analítica');
end
%%
%Gráfico semilog pw vs log t
if FILTRO9 == 3
 if FILTRO8 == 3
    plot7=plot(log(vector_t*tp/(vector_t+tp)),vector_Presion,'bo','LineWidth', 1.5);
    hold on
    plot8 = plot(log(vector_t*tp/(vector_t+tp)),Presion_analitica,'rs','LineWidth',1.5);

    grid on
    title( 'Prueba de Presion' , 'FontSize' ,12 ) 
    xlabel ( 'log (t [h]) ')
    ylabel ( 'Pw [lpc]')
    %legend('P.real','P.analítica');
 end
end
if FILTRO9 == 1
   if FILTRO8 == 3
    plot7=plot(log(vector_t*tp/(vector_t+tp)),vector_Presion,'bo','LineWidth', 1.5);
    hold on
    title( 'Prueba real' , 'FontSize' ,12 )
    grid on
    xlabel ( 'log (t [h]) ')
    ylabel ( 'Pw [lpc]')
   end
end
if FILTRO9 == 2
   if FILTRO8 == 3
    plot8 = plot(log(vector_t*tp/(vector_t+tp)),Presion_analitica,'rs','LineWidth',1.5);
    hold on
    title( 'Prueba analítica' , 'FontSize' ,12 )
    grid on
    xlabel ( 'log (t [h]) ')
    ylabel ( 'Pw [lpc]')
   end
end
%%
%Gráfico de la derivada
if FILTRO8==2%curva de la derivada
    if FILTRO2 == 2% Prueba PBU 
     
     deltap = curva_de_la_derivada_PBU(vector_presion,vector_t,longitud_vector_t,tp); 
     plot3=plot(log(vector_T),log(deltap),'bo','LineWidth', 1);
     hold on
     deltap = curva_de_la_derivada_PBU(Presion_analitica,vector_t,longitud_vector_t,tp); 
     plot4=plot(log(vector_T),log(deltap),'rs','LineWidth', 1);
     title( 'Curva de la derivada - Prueba PBU' , 'FontSize' ,12 )
     grid on
     xlabel ( '(delta t* tp)/(deltat + tp) [h] ')
     ylabel ( 'Delta P [lpc]')
     legend('P.real','P.analítica');
    end
    if FILTRO2 ==1% Prueba DD
      deltap = curva_de_la_derivada_DD(vector_presion,vector_t,longitud_vector_t); 
      plot5=plot(log(vector_T),log(deltap),'bo','LineWidth', 1);
      hold on
      deltap = curva_de_la_derivada_DD(Presion_analitica,vector_t,longitud_vector_t); 
      plot6=plot(log(vector_T),log(deltap),'rs','LineWidth', 1);
      title( 'Curva de la derivada - Prueba DD' , 'FontSize' ,12 )
      grid on
      xlabel ( 'Tiempo [h] ')
      ylabel ( 'Delta P [lpc]')
      legend('P.real','P.analítica');   
    end
end
%%
if FILTRO9 == 3
 if FILTRO8==4%Curva de la derivada & Rammey
    if FILTRO2 == 2% Prueba PBU 
     deltap = curva_de_la_derivada_PBU(Presion_analitica,vector_t,longitud_vector_t,tp); 
     plot4=plot(log((vector_T*tp)/(vector_T+tp)),log(deltap),'rs','LineWidth', 1);
     hold on
     deltap = Rammey(vector_presion,longitud_vector_t);
     plot7=plot(log((vector_t*tp)/(vector_t+tp)),log(deltap),'bo','LineWidth', 1);    
     hold on
     deltap = curva_de_la_derivada_PBU(vector_presion,vector_t,longitud_vector_t,tp); 
     plot3=plot(log((vector_T*tp)/(vector_T+tp)),log(deltap),'bo','LineWidth', 1);
     hold on 
     deltap = Rammey(Presion_analitica,longitud_vector_t); 
     plot8=plot(log((vector_t*tp)/(vector_t+tp)),log(deltap),'rs','LineWidth', 1);

     grid on
     title( 'Curva de la derivada  Prueba PBU' , 'FontSize' ,12 )
     xlabel ( '(delta t* tp)/(deltat + tp) [h] ')
     ylabel ( 'Delta P [lpc]')
    end
    
    if FILTRO2 ==1% Prueba DD
        
     deltap = curva_de_la_derivada_DD(Presion_analitica,vector_t,longitud_vector_t); 
     plot6=plot(log(vector_T),log(deltap),'rs','LineWidth', 1);
     hold on
     deltap = Rammey(Presion_analitica,longitud_vector_t); 
     plot8=plot(log(vector_t),log(deltap),'rs','LineWidth', 1);
     grid on 
     deltap = curva_de_la_derivada_DD(vector_presion,vector_t,longitud_vector_t); 
     plot5=plot(log(vector_T),log(deltap),'bo','LineWidth', 1);
     hold on
     deltap = Rammey(vector_presion,longitud_vector_t);
     plot7=plot(log(vector_t),log(deltap),'bo','LineWidth', 1);    
     hold on 
     
     title( 'Curva de la derivada - Prueba DD' , 'FontSize' ,12 )
     xlabel ( 'Tiempo [h] ')
     ylabel ( 'Delta P [lpc]')      
    end
 end
end
if FILTRO9 == 2% Resultados analíticos
  if FILTRO8==4%Curva de la derivada & Rammey
    if FILTRO2 == 2% Prueba PBU  
        deltap = curva_de_la_derivada_PBU(Presion_analitica,vector_t,longitud_vector_t,tp); 
        plot4=plot(log(vector_T),log(deltap),'rs','LineWidth', 1);
        hold on
        deltap = Rammey(Presion_analitica,longitud_vector_t);
        plot8= plot(log(vector_T),log(deltap),'rs','LineWidth', 1); 
     end
     title( 'Curva de la derivada - Prueba PBU' , 'FontSize' ,12 )
     %hold on
     grid on
     xlabel ( '(delta t* tp)/(deltat + tp) [h] ')
     ylabel ( 'Delta P [lpc]')
  end
end