function [B_T,no_T,ne_T,dno,dne] = fun_birefringence(material,lambda,T)
    if (strcmp(material,'Quartz'))
        lambda_um=lambda*10^(-3);
        %no
        no=sqrt(1.35362855744073+1.00314737858204.*lambda_um.^2./(lambda_um.^2-0.0106769269234773)+ ...,
        1.69460825115705.*(lambda_um.^2)./(lambda_um.^2-151.030393473906));
        lambda_ig_o_Quartz=1240/10.30*10^(-3);%E=1240/λ
        R_o_Quartz=lambda_um.^2./(lambda_um.^2-lambda_ig_o_Quartz.^2);
        dno_Quartz=(-61.1840*R_o_Quartz+43.999*R_o_Quartz.^2)./(2.*no).*10^(-6);
        no_T=no+dno_Quartz*(T-25);
        dno=dno_Quartz;
        %ne
        ne=sqrt(1.35388749649506+1.02903415330593.*lambda_um.^2./(lambda_um.^2-0.0108019305588061)+ ...,
        13.0793935135060.*(lambda_um.^2)./(lambda_um.^2-1209.93302155300)); 
        lambda_ig_e_Quartz=1240/10.30*10^(-3);%E=1240/λ
        R_e_Quartz=lambda_um.^2./(lambda_um.^2-lambda_ig_e_Quartz.^2);
        dne_Quartz=(-70.1182*R_e_Quartz+49.2875*R_e_Quartz.^2)./(2.*ne)*10^(-6);
        ne_T=ne+dne_Quartz*(T-25);
        dne=dne_Quartz;
        %birefingence
        B_T=ne_T-no_T;
    elseif strcmp(material,'Calcite')
        lambda_um=lambda*10^(-3);
        %no
        no=sqrt(1.73358749+0.96464345.*lambda_um.^2./(lambda_um.^2-0.0194325203)+ ...,
        1.82831454.*(lambda_um.^2)./(lambda_um.^2-120));
%         lambda_ig_o_calcite=1240/10.80*10^(-3);%E=1240/λ
        lambda_ig_o_calcite=0.152766988280853;
        R_o_calcite=lambda_um.^2./(lambda_um.^2-lambda_ig_o_calcite.^2);
%         dno_calcite=(-121.6890*R_o_calcite+122.4940*R_o_calcite.^2)./(2.*no).*10^(-6);
        dno_calcite=(-31.3565081593939*R_o_calcite+32.0188355454959	*R_o_calcite.^2)./(2.*no).*10^(-6);
        no_T=no+dno_calcite*(T-25);
        dno=dno_calcite;
        %ne
        ne=sqrt(1.35859695+0.82427830.*lambda_um.^2./(lambda_um.^2-0.0106689543)+ ...,
        0.14429128.*(lambda_um.^2)./(lambda_um.^2-120)); 
%         lambda_ig_e_calcite=1240/9.05*10^(-3);%E=1240/λ
        lambda_ig_e_calcite=0.127015496568904;
        R_e_calcite=lambda_um.^2./(lambda_um.^2-lambda_ig_e_calcite.^2);
%         dne_calcite=(12.7011*R_e_calcite+20.4803*R_e_calcite.^2)./(2.*ne)*10^(-6);
        dne_calcite=(-5.81366738478945*R_e_calcite+35.6342668226812*R_e_calcite.^2)./(2.*ne)*10^(-6);
        ne_T=ne+dne_calcite*(T-25);
        dne=dne_calcite;
        %birefingence
        B_T=ne_T-no_T;
    elseif strcmp(material,'MgF2')
        lambda_um=lambda*10^(-3);
        %no
        no=sqrt(1.30618931+0.58007346.*lambda_um.^2./(lambda_um.^2-0.0077355208)+ ...,
        2.24075551.*(lambda_um.^2)./(lambda_um.^2-550));
        lambda_ig_o_MgF2=1240/13.10*10^(-3);%E=1240/λ
        R_o_MgF2=lambda_um.^2./(lambda_um.^2-lambda_ig_o_MgF2.^2);
        dno_MgF2=(-37.2043*R_o_MgF2+39.3186*R_o_MgF2.^2)./(2.*no).*10^(-6);
        no_T=no+dno_MgF2*(T-25);
        dno=dno_MgF2;
        %ne
        ne=sqrt(1.29108974+0.62728316.*lambda_um.^2./(lambda_um.^2-0.0075322929)+ ...,
        2.41840547.*(lambda_um.^2)./(lambda_um.^2-550)); 
        lambda_ig_e_MgF2=1240/15.50*10^(-3);%E=1240/λ
        R_e_MgF2=lambda_um.^2./(lambda_um.^2-lambda_ig_e_MgF2.^2);
        dne_MgF2=(-56.7859*R_e_MgF2+57.3986*R_e_MgF2.^2)./(2.*ne)*10^(-6);
        ne_T=ne+dne_MgF2*(T-25);
        dne=dne_MgF2;
        %birefingence
        B_T=ne_T-no_T;
    else
        disp('请输入Quartz、Calcite、MgF2三种材料中的一种');
    end
end