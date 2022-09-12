close all;
clearvars -except direction* v_result_average* displacement_* msd*;
run=8;

Estrength_initialize = [0.95];
noise_initialize = [0.05];
ra_initialize = [2];
cell_number=[25 35 45 65];
Electric_Noise = struct();
Electric_Noise.direction(1:length(Estrength_initialize), 1:length(noise_initialize)) = zeros(length(Estrength_initialize), length(noise_initialize));

for Estrength_loop=1:length(Estrength_initialize)
    for cellnum=1:length(cell_number)
        for run_iteration=1:run
            close all; 
            clearvars -except direction* v_result_average* ang_result tot_v run run_iteration displacement_* ra_* Estrength_* noise_* Electric_Noise msd* cell_number* cellnum;
            lbox=150.;
            tot_time = 130;
            ncell = cell_number(cellnum);
            t_relax = 15;
            vs_init=0.15; 
            ve_init = 0.004;  
            vrep_init=0.3;
            E_strength = Estrength_initialize(Estrength_loop);
            eta=0.05; 
            eta_efield =0;
            k=0.3;
            dt=1; 
            scale_efield = 2;
            mag_E=0.5;
      
            r=ra_initialize;
            r_rep= 2;
            r_efield = 10;
            mindex=0; 
            
            filename_movie = sprintf('traj_estrength%1.2f_e%.3f_vo%.1f_vrep%.2f_eta%.2f_alignr%d_iteration%d', E_strength,ve_init, vs_init, vrep_init, eta,r,run_iteration);
            filename_traj = sprintf('traj_estrength%1.2f_e%.3f_vo%.1f_vrep%.2f_eta%.2f_alignr%d_iteration%d', E_strength,ve_init, vs_init, vrep_init, eta,r,run_iteration);
            
            celll=[1:ncell]; 
            onebd=ones(1,ncell)'; 
            %cell arrays position and angle 
            %plot cell trace
            xy_cell_i = nan(tot_time,2);
            particles=struct();
            distance = struct();
            
            cell=figure; 
            axis([0 lbox 0 lbox])
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',12);
            axis('square') 
            hold on 

            init_R = lbox/8;
            init_theta = 2*pi*rand(ncell,1);
            init_r = init_R*sqrt(rand(ncell,1));
            xb = lbox/2 + init_r.*cos(init_theta);
            yb = lbox/2 + init_r.*sin(init_theta);

            ang=pi.*2.*(rand(ncell,1));
             
            vxb=vs_init.*cos(ang); 
            vyb=vs_init.*sin(ang); 
             
            %outer loop 
       
            %Discretize spatial grid for E-field
            [x,y] = meshgrid([0.1:2:lbox+0.5]);
            x_efield = reshape(x,length(x)^2,1);
            y_efield = reshape(y,length(y)^2,1);
            x_efield(length(x_efield)+1:length(onebd))=0;
            y_efield(length(y_efield)+1:length(onebd))=0;
            u=(E_strength)*x./(x);
            v= E_strength*1e-9*y;
          
            u_efield = reshape(u,length(u)^2,1);
            v_efield = reshape(v,length(v)^2,1);
            u_efield(length(u_efield)+1:length(onebd))=0;
            v_efield(length(v_efield)+1:length(onebd))=0;
            
            %for plotting efields 
            x_efield_plot = x_efield(1:1:end);
            y_efield_plot = y_efield(1:1:end);
            u_efield_plot = u_efield(1:1:end);
            v_efield_plot = v_efield(1:1:end);
            %---------
            tot_uv_efield = u_efield+v_efield;
 
            ang_efield = atan(v_efield./u_efield);
            ang_efield(isnan(ang_efield))=0;
            ang_efield(length(ang_efield)+1:length(onebd)) = 0;
            
            onebd_efield=ones(1,length(x_efield))';
            
            e_freq = linspace(0,2*pi, tot_time);
            Tot_energy= zeros(tot_time,1);
            phi = zeros(tot_time,1);
            
            curve = animatedline;
            for nsteps=1:tot_time; 
                nsteps 
                if(nsteps<t_relax)
                    vs = 0; ve=0; vrep = vrep_init;
                   
                else ve=ve_init; vrep = vrep_init; vs= vs_init;%*rand(ncell,1);
                end 
                
                xb=xb+vxb.*dt; 
                yb=yb+vyb.*dt; 
                 
                particles.x(nsteps,:) = xb';
                particles.y(nsteps,:) = yb';
                
                for cell1=1:ncell
                    %periodic 
                    if(xb(cell1)<0);xb(cell1)=xb(cell1)+lbox; end 
                    if (yb(cell1)<0);yb(cell1)=yb(cell1)+lbox;end 
                    if (xb(cell1)>lbox);xb(cell1)=xb(cell1)-lbox;end 
                    if (yb(cell1)>lbox);yb(cell1)=yb(cell1)-lbox;end 
                     
                    %find mean angle of neigbours (include cell1) 
                    sep(1:ncell)=sqrt((onebd.*xb(cell1)-xb(1:ncell)).^2+... 
                        (onebd.*yb(cell1)-yb(1:ncell)).^2); 
                    sep_efield(1:length(x_efield))=(onebd_efield.*xb(cell1)-x_efield).^2+... 
                        (onebd_efield.*yb(cell1)-y_efield).^2; 
                    
                    dist_x(1:ncell)=(-onebd.*xb(cell1)+ xb(1:ncell));
                    dist_y(1:ncell)=(-onebd.*yb(cell1)+yb(1:ncell));
                    ang_force=atan2(dist_y(find(sep<r_rep & sep>0)), dist_x(find(sep<r_rep & sep>0) ));
                    ang_force_mean(cell1)=mean(ang_force);
                    %seperate(1:ncell)=sqrt(dist_x(1:ncell).^2+dist_y(1:ncell).^2);
                    %forcex(cell1)=mean(-k.*(1./(1+(sep(sep<r_rep & sep>0)')./r_rep)).*cos(ang_force'));
                    %forcey(cell1)=mean(-k.*(1./(1+(sep(sep<r_rep & sep>0)')./r_rep)).*sin(ang_force'));
                    forcex(cell1)=mean(-k.*((2-(sep(sep<r_rep & sep>0)')./r_rep)).*cos(ang_force'));
                    forcey(cell1)=mean(-k.*((2-(sep(sep<r_rep & sep>0)')./r_rep)).*sin(ang_force'));
                    %avforcex=-k.*(min_sep)*cos(ang_force);
                    %avforcey=-k.*(min_sep)*sin(ang_force);
                    forcex(isnan(forcex))=0;
                    forcey(isnan(forcey))=0;
                    mag_force(cell1)=sqrt(forcex(cell1).^2+forcey(cell1).^2);
                    ang_force2(cell1)=atan2(forcey(cell1),forcex(cell1) );
                    
                    nearang=ang(sep<r);  %'r' this is alignment radius 
                    mang(cell1)=mean(nearang);
                    
                    nearang_efield = atan2((v_efield(sep_efield<r_efield)), u_efield(sep_efield<r_efield)) ;
                    nearmag_efield_x = u_efield(sep_efield<r_efield);
                    nearmag_efield_y = v_efield(sep_efield<r_efield);
                    mnearmag_efield_x(cell1) = mean(nearmag_efield_x);
                    mnearmag_efield_y(cell1) = mean(nearmag_efield_y);
                    mang_efield(cell1) = mean(nearang_efield);
                    mang_cell_efield(cell1) = (1*mang(cell1)+mang_efield(cell1))/2;
                end 
                
                ang_force2(isnan(ang_force2))=0;
                ang_force_mean(isnan(ang_force_mean))=0;
                ang_force2_trans=ang_force2';
    
                ang=mang' + eta*(rand(ncell,1)-0.5)*2*pi + 0*mang_cell_efield'+0*ang_force2' ;
                ang_efield = mang_efield' + eta_efield*(rand(ncell,1)-0.5)*2*pi;
                
                near_efield_x = abs(mnearmag_efield_x);
                near_efield_y = abs(mnearmag_efield_y);
    
                
                 vxb=1*vs.*cos((ang+0*ang_efield)./1)+ (vrep.*mag_force.*(cos(ang_force2)))' - ve*near_efield_x*cos(ang_efield) ;
                 vyb=1*vs.*sin((ang+0*ang_efield)./1)+ (vrep.*mag_force.*(sin(ang_force2)))' - ve*near_efield_y*sin(ang_efield);
               
                vx_norm = vxb./norm(vxb);
                vy_norm = vyb./norm(vyb);
                v_result = [vxb vyb];
                v_result_norm = sqrt(diag(v_result * v_result'));
    
                Tot_energy(nsteps,1) = (sum(forcex(:))^2+sum(forcey(:))^2)/ncell;
                
                xy_cell_i(nsteps,1:2) = [xb(1) yb(1)]; 
    
                cla 
                set(gcf,'doublebuffer','on')  
                pcol = pcolor(x,y,1*ones(size(u)));
                pcol.FaceColor = 'interp';
                pcol.EdgeColor = 'interp';
                colorMap = [0.9290, 0.7940, 0.7250];
                colormap(colorMap);
                
                hold on;
                skip_nth =14;
                cell.WindowState = 'maximized';
            end %end of time loop
            
            %MSD calculation
            for time=1:tot_time
                distance.d(time,:)=sqrt(particles.x(time,:).*particles.x(time,:)-particles.y(time,:).*particles.y(time,:));
                for cell1=1:ncell
                    eval(['msd' strrep(num2str(E_strength), '.', 'p') '_e' strrep(num2str(ve_init), '.', 'p') '_vo' strrep(num2str(vs_init), '.', 'p') '_vrep' strrep(num2str(vrep_init), '.', 'p') '_eta' strrep(num2str(eta), '.', 'p') '_alignr' strrep(num2str(num2str(r,'%.5f')), '.', 'p') '(time,1)'  '=(sum((distance.d(time,cell1)-distance.d(1,cell1))*(distance.d(time,cell1)-distance.d(1,cell1))))/(ncell);']);
            
                end
            end
            
            
            
            for t=1:ncell
                ang_result(t,1*(run_iteration-1)+1:1*(run_iteration-1)+1) = atan2d(particles.y(end,t)-particles.y(1,t),particles.x(end,t)-particles.x(1,t)) + 360*((particles.y(end,t)-particles.y(1,t))<0);
                tot_v(t,1*(run_iteration-1)+1:1*(run_iteration-1)+1) = sqrt((particles.y(end,t)-particles.y(1,t)).^2+(particles.x(end,t)-particles.x(1,t)).^2)/(tot_time-t_relax);
            end
            
            close all;
            eval(['direction_estrength' strrep(num2str(E_strength), '.', 'p') '_e' strrep(num2str(ve_init), '.', 'p') '_vo' strrep(num2str(vs_init), '.', 'p') '_vrep' strrep(num2str(vrep_init), '.', 'p') '_eta' strrep(num2str(eta), '.', 'p') '_alignr' strrep(num2str(num2str(ncell,'%.5f')), '.', 'p') '(1,run_iteration)'  '=mean(cos((ang_result(:,run_iteration)).*(pi/180)));']);
            eval(['v_result_average' strrep(num2str(E_strength), '.', 'p') '_e' strrep(num2str(ve_init), '.', 'p') '_vo' strrep(num2str(vs_init), '.', 'p') '_vrep' strrep(num2str(vrep_init), '.', 'p') '_eta' strrep(num2str(eta), '.', 'p') '_alignr' strrep(num2str(num2str(r,'%.5f')), '.', 'p') '(1,run_iteration)'  '=mean(((tot_v(:,run_iteration))));']);
           
            end %ends the iteration loop
end %ends the 'r_a' loop
end
Electric_Noise.direction(1:end, 1:end)=Electric_Noise.direction(1:end, 1:end)/run;
ang_result = reshape(ang_result, ncell*run_iteration,1);
tot_v = reshape(tot_v,ncell*run_iteration,1);
 
 for en=1:(tot_time/4) % given that tot_time is a multpile of 20 or your choice of no. of rows to average.
    a=(en-1)*4+1;
    b=(en-1)*4+4;
    meandata(en,1) = mean(Tot_energy(a:b));
    meanphi(en,2) = mean(phi(a:b));
    meanphi(en,1) = (a+b)/2;
 end
 
