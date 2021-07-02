%dim eslabones
R=2.102;V=1.4;B=1.2;W=0.327; r=0.15;%largo eslabones 
a1=0.85;a2=0.7;a3=0.4; %espesor esl

%masas e inercias
m1=1333;m2=1067;m3=267;m4=310;%c1=1;c2=1;c3=1;c4=1;
dcm34=((m3*W/2)+(m4*B))/(m3+m4);

J1y=(m1*(R^2+a1^2))/12;J2y=(m2*(V^2+a2^2))/12;
J3y=(m1*(W^2+a3^2))/12; J4y=(m4*(3*(r^2)+(2*(B-W)^2)))/12;
J34y=J3y+(m3*(dcm34-(W/2))^2)+J4y+(m4*(B-(dcm34))^2);
J4z=(m4*r^2)/2; %arreglar J34y
g=9.8;


%valores de K
k1=3.8*(10^6);k2=6.6*(10^6);k3=3.9*(10^6);k4=5.6*(10^5); %Nm/rad Ktheta
%M K C equivalentes
Meq=[J1y+((m1*(R^2))/4)+m2*(R^2) 0 (m3+m4)*(R*dcm34) (m3+m4)*(R^2);
    0 J2y+(m2*(V^2)/2)+((V^2)*(m3+m4)) -(m3+m4)*V*dcm34 0;
    (m3+m4)*(R*dcm34) -(m3+m4)*V*dcm34 J34y+((m3+m4)*(dcm34^2)) 0;
    0 0 0 J4z];
Keq=[(k1+k2)-(m1*g*R/2)-(m2*g*R)-((m3+m4)*g*(R-dcm34)) -k2 0 0;
    -k2 (k2+k3) -k3 0;
    0 -k3 ((m3+m4)*g*dcm34)+k3 0;
    0 0 0 k4];
%Ceq=[c1+c2 -c2 0 0;-c2 c2+c3 -c3 0;0 -c3 c3 0;0 0 0 c4];

%xita
%xita=0.02;
xita=[0.02;0.03;0.04;0.043];
%valores de contorno
x0=[0;0;0;0];%posiciones iniciales
%v0=[1.0123;1.0123;1.0123;1.0123];%velocidades iniciales, distinto a 0
v0=[0;1.0123;0;0];
gl=4;      %grados de libertad
%tiempo sim
dt=0.02;
tf=10;
t=0:dt:tf; 

%autovectores-autovalores
[autovec,autoval]=eig(Keq,Meq);
w=diag(autoval).^0.5;   %vector de frecuencias naturales
for i=1:gl
wd(i)=w(i)*(1-xita(i)^2)^0.5; 
end
I=eye(size(Meq));
M=autovec'*Meq*autovec;  % masa modal
q0=M\I*autovec'*Meq*x0; % desplazamientos iniciales en coords modales
qprima0=M\I*autovec'*Meq*v0; %velocidades iniciales en coords modales 

%descomposicion Modal
N=length(t);
q=[0];
t=t';
for i=1:gl
    for j=1:N
    q(i,j)=exp(-xita(i)*w(i)*t(j))*((cos(wd(i)*t(j))+(xita(i)/(1-xita(i)^2)^0.5)*sin(wd(i)*t(j))))*(q0(i))+(1/(wd(i)))*exp(-xita(i)*w(i)*t(j))*sin(wd(i)*t(j))*(qprima0(i));
    end
end
%tercer termino Duhamel regla de los trapecios
format long;
mm=diag(Meq);
F=[0;0;0;0];%F=zeros(length(t),length(t),length(t),length(t));
Q=autovec'*F;
A=zeros(gl,length(t));
B=A;
for k=1:gl
    for i=2:N
        j=i-1;
        yc(k,j)=0*cos(wd(k)*t(j)); %Q(k,j)            
        ys(k,j)=0*sin(wd(k)*t(j));
        A(k,i)=A(k,j)*exp(-xita(k)*dt*w(k))+dt/(2*wd(k))*(yc(k,j)*exp(-xita(k)*dt*w(k))+yc(k,j));
        B(k,i)=B(k,j)*exp(-xita(k)*dt*w(k))+dt/(2*wd(k))*(ys(k,j)*exp(-xita(k)*dt*w(k))+ys(k,j));
    end
end
for k=1:gl
    for j=1:N
        qInt(k,j)=mm(k)*(A(k,j)*sin(wd(k)*t(j))-B(k,j)*cos(wd(k)*t(j)));
    end
end
%q=q+qInt;
x=autovec*q;

%% Desplazamiento para cada Articulacion
figure('Name','Posiciones');
subplot(4,1,1);plot(t,x(1,:),'-r');xlabel("tiempo");ylabel("Desp. [rad]");title("1er GDL");legend("Descomp. modal");
subplot(4,1,2);plot(t,x(2,:),'-r');xlabel("tiempo");ylabel("Desp. [rad]");title("2do GDL");legend("Descomp. modal");
subplot(4,1,3);plot(t,x(3,:),'-r');xlabel("tiempo");ylabel("Desp. [rad]");title("3er GDL");legend("Descomp. modal");
subplot(4,1,4);plot(t,x(4,:),'-r');xlabel("tiempo");ylabel("Desp. [rad]");title("4to GDL");legend("Descomp. modal");

%velocidad
v1=diff(x(1,:))/dt;
v2=diff(x(2,:))/dt;
v3=diff(x(3,:))/dt;
v4=diff(x(4,:))/dt;
v=[v1;v2;v3;v4];
aux=length(t)-1;
ti=(1:aux)*dt;
figure('Name','Velocidades');
subplot(4,1,1);plot(ti,v(1,:),'-m');xlabel("tiempo");ylabel("Vel [rad/s]");title("1er GDL");legend("Descomp. modal");
subplot(4,1,2);plot(ti,v(2,:),'-m');xlabel("tiempo");ylabel("Vel [rad/s]");title("2do GDL");legend("Descomp. modal");
subplot(4,1,3);plot(ti,v(3,:),'-m');xlabel("tiempo");ylabel("Vel [rad/s]");title("3er GDL");legend("Descomp. modal");
subplot(4,1,4);plot(ti,v(4,:),'-m');xlabel("tiempo");ylabel("Vel [rad/s]");title("4to GDL");legend("Descomp. modal");
%aceleraciones
aux=aux-1;
ti=(1:aux)*dt;
a1=diff(v(1,:))/dt;
a2=diff(v(2,:))/dt;
a3=diff(v(3,:))/dt;
a4=diff(v(3,:))/dt;
a=[a1;a2;a3;a4];
figure('Name','Aceleraciones');
subplot(4,1,1);plot(ti,a(1,:),'-b');xlabel("tiempo");ylabel("Ac. [rad/s^2]");title("1er GDL");legend("Descomp. modal");
subplot(4,1,2);plot(ti,a(2,:),'-b');xlabel("tiempo");ylabel("Ac. [rad/s^2]");title("2do GDL");legend("Descomp. modal");
subplot(4,1,3);plot(ti,a(3,:),'-b');xlabel("tiempo");ylabel("Ac. [rad/s^2]");title("3er GDL");legend("Descomp. modal");
subplot(4,1,4);plot(ti,a(4,:),'-b');xlabel("tiempo");ylabel("Ac. [rad/s^2]");title("4to GDL");legend("Descomp. modal");


%calculo matriz c
for i=1:gl
    for j=1:gl
        if(j==i)
           cm(i,j)=2*w(j)*xita(j);
        end
    end
end

c=Meq*autovec'*cm*autovec*Meq;
%c=autovec'*cm*autovec;
A=[1 1 0 0,
    0 1 1 0,
    0 0 1 0
    0 0 0 1];
B=[c(1,1), c(2,2), c(3,3), c(4,4)];
lsqr(A,B(:))
%mldivide(Ceq,c)