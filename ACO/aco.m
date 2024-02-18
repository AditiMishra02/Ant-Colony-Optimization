clear all
close all
clc

%% GPU config
% gpu = gpuDevice;
% gpu(1);

tic
%% Main configuration values for this simulation

dataset.nodeNo = 9; %Number of nodes
ACOnodeNo = dataset.nodeNo;
dataset.nodePosition(1,:) = [1 50 50]; %(Sender node fixed position)
dataset.nodePosition(2,:) = [2 900 900]; %(Receiver node fixed position)
dataset.NeighborsNo = 5;
dataset.range = 500; %Tolerance distance to became neighbor of one node (Euclidean distance based)
dataset.atenuationFactor = 1.8; %Atenuation factor in freespace - ranges from 1.8 to 4 due environment
dataset.minEnergy = 80; % Mw - Miliwatts (70% energy)
dataset.maxEnergy = 100; % Mw - Miliwatts (Full energy (100%) - 1 mAh charge capacity within 1 Volt energy)
dataset.energyconsumptionperCicle = 0.85;
dataset.energyrecoveryperCicle = 0.2;
dataset.minenergyfactor = 0.18;
dataset.maxenergyfactor = 0.2;
STenergy=inf; 
packet=0;
iterationcounter=1;
plotgraphs=1; %Choose 1 for "yes" or 0 for "no" if you want to plot graphs or no (Better performance if no)
reprodutibily = 0; %1 = yes (always generate same random numbers) (0) for no reprodutibility (Different random numbers every code execution);


% Node position sortition
if reprodutibily == 0
    rng('shuffle');
else
    rng('default');
end
for a = 3 : dataset.nodeNo
    
   dataset.nodeId = a; 
   garbage.x = randi([1 900]); %Xpos sortition
   garbage.y = randi([1 900]); %Ypos sortition
   dataset.nodePosition(a,:) = [dataset.nodeId garbage.x garbage.y]; %NodeID, X and Y position into nodePosition table
   
end

% Euclidean Distance calc from one node to all others

for i = 1 : dataset.nodeNo
    for j = 1: dataset.nodeNo
        garbage.x1 = dataset.nodePosition(i,2); 
        garbage.x2 = dataset.nodePosition(j,2); 
        garbage.y1 = dataset.nodePosition(i,3); 
        garbage.y2 = dataset.nodePosition(j,3);
        
        dataset.euclidiana(i,j) = sqrt(  (garbage.x1 - garbage.x2) ^2 + (garbage.y1 - garbage.y2)^2  ); 
        
    end
end

% Edges matrix definition due "range" variable value

dataset.weights = lt(dataset.euclidiana,dataset.range);

% Graph construction

G=graph(dataset.weights,'omitselfloops'); %Graph creation based on adjacency matrix (Edges matrix) built above

% Euclidean distance extraction for all existente end-to-end formed by
% "distance tolerance" (range variable value)

for a = 1 : height(G.Edges)
    garbage.s = G.Edges.EndNodes(a,1);
    garbage.t = G.Edges.EndNodes(a,2);
    garbage.Z(a,:) = dataset.euclidiana(garbage.s,garbage.t);
 end
G.Edges.Euclidiana = garbage.Z(:,1);

%Initial energy sortition (from 70% to 100% - minEnergy and maxEnergy variable valeu) 

[dataset.nodePosition(:,4)] = dataset.maxEnergy -(dataset.maxEnergy-dataset.minEnergy)*rand(dataset.nodeNo,1);
dataset.nodePosition(1:2,4)=STenergy;

%All "G" (Graph object) based nodes degree to use as "node processing
%status overload" (more connections, busier!)

for a = 1: length(dataset.nodePosition(:,1))
   
    dataset.nodePosition(a,5) = degree(G,dataset.nodePosition(a,1));
    
end

% Pathloss calc of each Edges based in a freespace (1.8 factor)

[G.Edges.Pathloss] = (10*dataset.atenuationFactor)*log10(G.Edges.Euclidiana);

%End points coordinates and extra data migration to G object

for a = 1 : height(G.Edges)
	garbage.Sourcenode = G.Edges.EndNodes(a,1);
	garbage.Targetnode = G.Edges.EndNodes(a,2);
	G.Edges.SourcenodeXpos(a) = dataset.nodePosition(garbage.Sourcenode,2);
	G.Edges.SourcenodeYpos(a) = dataset.nodePosition(garbage.Sourcenode,3);
	G.Edges.TargetnodeXpos(a) = dataset.nodePosition(garbage.Targetnode,2);
	G.Edges.TargetnodeYpos(a) = dataset.nodePosition(garbage.Targetnode,3);
    G.Edges.ActiveEdge(a) = 1;
end

%G Edges cumulative energy calc
for a= 1: height(G.Edges)
	
	G.Edges.Edgeenergy(a,:)=dataset.nodePosition(G.Edges.EndNodes(a,1),4)+dataset.nodePosition(G.Edges.EndNodes(a,2),4);

end

%Cumulative edges energy matrix creation
dataset.edgesenergymatrix=zeros(dataset.nodeNo);
for a = 1 : height(G.Edges)
	dataset.edgesenergymatrix(G.Edges.EndNodes(a,1),G.Edges.EndNodes(a,2))=1/G.Edges.Edgeenergy(a);
	dataset.edgesenergymatrix(G.Edges.EndNodes(a,2),G.Edges.EndNodes(a,1))=1/G.Edges.Edgeenergy(a);
end
ACOedgesenergymatrix=dataset.edgesenergymatrix;

% Graph objects plot

if plotgraphs == 1
    figure('units','normalized','innerposition',[0 0 1 1],'MenuBar','none')
    subplot(1,1,1) %1,3,1 (Line number,collumn number, graph id) - if you want to show more than 1 graph in same windows
    garbage.Xmax = 1500;
    garbage.Xmin = 0;
    garbage.Ymax = 1500;
    garbage.Ymin = 0;
    p = plot(G,'XData',(dataset.nodePosition(:,2)),'YData',(dataset.nodePosition(:,3))); 
    line(dataset.nodePosition(1:2,2),dataset.nodePosition(1:2,3),'color','green','marker','o','linestyle','none','markersize',50)
    garbage.ax = gca;
    garbage.ax.XAxis.TickValues = 0:100:1000;
    garbage.ax.YAxis.TickValues = 0:100:1000;
    grid on
    hold on
    title(['Original WSN | ','Nodes number: ',num2str(dataset.nodeNo),' | Nodes range: ', num2str(dataset.range)])
    figname=[strcat('FIG-Original-NO-ACO-simulation-',num2str(dataset.nodeNo),'nodes-',datestr(now,'mm-dd-yyyy-HH-MM'),'.png')];
    saveas(gcf,figname)
    pause(2)
end

% Rebuild (clone) original graph for ACO workload

Gaco = G; %ACO Graph creation cloned from original G object

%Creation of some stuff needed
garbage.deadnodelist=[];
garbage.deadnodeneighbors=[];

%% finding routes

% Finding shortest path route
G2 = shortestpathtree(G,1,2);

%Clone dataset before changes to use in ACO script
ACOdataset = dataset;
pause(2)

%Prepare output filename to save the reports
name=[strcat('NO-ACO-simulation-',num2str(dataset.nodeNo),'nodes-',datestr(now,'mm-dd-yyyy-HH-MM'),'.txt')];
fileID = fopen(fullfile(pwd,name),'w');
fprintf(fileID,'%6s %20s %20s %20s %20s\r\n','|NodeNo|','|No ACO Scene|','|Hops|','|Packets sent|','|Dead node|');


%% Initialize path existance test for loops


while ~isempty(G2.Edges)
    G2 = shortestpathtree(G,1,2);
    iterationcounter=iterationcounter+1;
    % Test if there is connection between node 1 and 2. If not, terminate!
    if isempty(G2.Edges)
        break
    end
    
    %Find edges found by shorthestpathtree in main G object
    for a = 1 : height(G2.Edges)
        source = G2.Edges.EndNodes(a,1);
        target = G2.Edges.EndNodes(a,2);
        garbage.edgesrow(a,:) = findedge(G,source,target);
    end

    %Find nodes involved in routing event
    garbage.routingnodes = unique(G2.Edges.EndNodes);
    garbage.routepath = shortestpath(G,1,2);
    %shortestpath(G,1,2)
            
    %Construct localization dataset to nodes involved in routing event for plot effects
    for a = 1 : length(garbage.routingnodes)
        garbage.b=garbage.routingnodes(a,1);
        dataset.routingnodesPosition(a,:)=dataset.nodePosition(garbage.b,:);
    end
        
    %Find nodes not involved in routing and its energy to energy increase
    [garbage.Outroutenodes, garbage.Outroutenodesidx] = setdiff(dataset.nodePosition(:,1),garbage.routingnodes(:,1),'stable');
    for a = 1 : length(garbage.Outroutenodes)
        garbage.b=garbage.Outroutenodes(a,1);
        garbage.Outroutenodes(a,2)=dataset.nodePosition(garbage.b,4);
    end
    
    
    %Code block for energy usage (decrease) and packet send count
    
    while min(dataset.nodePosition(:,4))>0 
        randomenergyfactor=(dataset.maxenergyfactor-dataset.minenergyfactor).*rand(1,1);

        for a = 1 : length(garbage.routingnodes)
            node=garbage.routingnodes(a,1);
            dataset.nodePosition(node,4)=dataset.nodePosition(node,4)-dataset.energyconsumptionperCicle^randomenergyfactor+dataset.energyrecoveryperCicle^randomenergyfactor;
            
        end
        packet=packet+1;
        if min(garbage.Outroutenodes(:,2))<90
                [garbage.Lowenergyidx] = find(garbage.Outroutenodes(:,2)<=90);
                for b = 1 : length(garbage.Lowenergyidx)
                    node=garbage.Outroutenodes(garbage.Lowenergyidx(b,1),1);
                    garbage.Outroutenodes(b,2)=garbage.Outroutenodes(b,2)+dataset.energyrecoveryperCicle^randomenergyfactor;
                    dataset.nodePosition(node,4)=garbage.Outroutenodes(b,2);
                end
        end
    end
    
    %Find dead node ID to discorver its neighbors and disable relative edges 
    [garbage.deadnoderow] = find(dataset.nodePosition(:,4)<=0);
    for a = 1 : length(garbage.deadnoderow)
        deadnode=garbage.deadnoderow(a,1);
        for b = 1 : height(G.Edges)
            if ismember(G.Edges.EndNodes(b,1),deadnode) == 1 || ismember(G.Edges.EndNodes(b,2),deadnode) == 1
                G.Edges.ActiveEdge(b)=0;
                deadnode;
                pause(2)
            end
        end
    end
    garbage.deadnodelist(length(garbage.deadnodelist)+1,1)=deadnode;
    [garbage.deadedgerow]=find(G.Edges.ActiveEdge==0);

    %Mark node energy as NaN
    [garbage.deaddatasetnoderow]=find(dataset.nodePosition(:,4)<=0);
    for a = 1 : length(garbage.deaddatasetnoderow)
        b=garbage.deaddatasetnoderow(a,1);
        dataset.nodePosition(b,4)=NaN;
    end
    
    %Get number of hops between source and destination for this routing
    %session
    hopsnumber=length(garbage.routingnodes);

    %Informative section in command line - If need to disable plot loop to increase
    %processing speed
    
    msg=['Without ACO - Scene #: ',num2str(iterationcounter),' | Hops : ',num2str(hopsnumber),' | Packets sent: ', num2str(packet),' | Dead node: ', num2str(deadnode),' | Routing nodes: ', num2str(garbage.routepath)];
    disp(msg)
    
    if plotgraphs == 1
    %plot for every dead edges iteration's result
    figure('units','normalized','innerposition',[0 0 1 1],'MenuBar','none')
    p = plot(G,'XData',(dataset.nodePosition(:,2)),'YData',(dataset.nodePosition(:,3))); 
    line(dataset.nodePosition(1:2,2),dataset.nodePosition(1:2,3),'color','green','marker','o','linestyle','none','markersize',50)
    garbage.ax = gca;
    garbage.ax.XAxis.TickValues = 0:100:1000;
    garbage.ax.YAxis.TickValues = 0:100:1000;
    hold on
    
    %Plot all dead nodes in red
    for a = 1 : length(garbage.deadnodelist)
        garbage.b=garbage.deadnodelist(a,1);
        scatter(dataset.nodePosition(garbage.b,2),dataset.nodePosition(garbage.b,3),'MarkerFaceColor','red');
    end
    
    %Title for every iteration's result
    title(['WSN shortest path: ',num2str(iterationcounter),' |hops: ',num2str(hopsnumber),' |Packets sent: ',num2str(packet),' |Dead node: ',num2str(deadnode),' |Router nodes: ', num2str(garbage.routepath)])
    grid on
    pause(2)
       
    %Mark nodes involved in route with green color
    scatter(dataset.routingnodesPosition(:,2),dataset.routingnodesPosition(:,3),'MarkerFaceColor','green');
    
    figname=[strcat('FIG-NO-ACO-simulation-',num2str(dataset.nodeNo),'nodes-','scene-',num2str(iterationcounter),'-date-',datestr(now,'mm-dd-yyyy-HH-MM'),'.png')];
    saveas(gcf,figname)
    
    pause(0.2) 
    
    
    end
    
    
    %Remove dead edges from graph
    G = rmedge(G,garbage.deadedgerow(:,1));


    %Colect data to append to a external file for later usage - file open
    %command line in the begining

    report1=[dataset.nodeNo;iterationcounter;hopsnumber;packet;deadnode];
    fprintf(fileID,'%6.0f %20.0f %20.0f %20.0f %20.0f\r\n',report1);
    
    %Clear router nodes position to avoid erroneous plot
    clear dataset.routingnodesPosition
end

fclose(fileID); %Close external data collector file

 if plotgraphs == 1
%Plot all dead nodes in red
    for a = 1 : length(garbage.deadnodelist)
        garbage.b=garbage.deadnodelist(a,1);
        scatter(dataset.nodePosition(garbage.b,2),dataset.nodePosition(garbage.b,3),'MarkerFaceColor','red');
    end
 end

%Title for the last plot
try
    title(['WSN shortest path: ',num2str(iterationcounter),' |hops: ',num2str(hopsnumber),' |Packets sent: ',num2str(packet),' |Dead node: ',num2str(deadnode),' |Router nodes: ', num2str(garbage.routepath),'{\color{red} - ALL ROUTES UNAVAILABLE}'])
catch
    disp('NO ROUTES BETWEEN SOURCE (NODE1) AND TARGET (NODE2)')
end

%Message if there is no path between source and target even in the first iteration
disp('NO ROUTES BETWEEN SOURCE (NODE1) AND TARGET (NODE2)')

%Call aco.m file to run same scenario within aco script
toc
run aco.m



disp('################ STARTING ACO SCRIPT ################')

%% Main configuration values for ACO simulation

ACOsourcenode = 1;
ACOtargetnode = 2;
ACOmaxAntiteration = 100;
ACOantsNo = 21;
ACOpheromone = 1; %tau - pheromone concentration - most papers use 0 as its value
% ACOpathloss = 0.2; %eta - pathloss of an edge for aco calc
% ACOedgeenergy = 0.5; %Energy weight for aco calc
ACOevapcoeficient = 0.2; % 5 percent - (rho) pheromone evaporation rate
ACOalpha = -0.1; %Pheromone exponent
ACObeta = 0.1; %Energy exponent
ACOomega = 1; %Pathloss expoent
ACOpacket=0;
ACOiterationcounter=1;
ACObestfitness=inf;
ACOenergycicles=1;
ACObesttourhistory=[];


%Replot original graph that will be used by ACO - just for comparison
%Gaco=graph(ACOedgesenergymatrix,'omitselfloops');
if plotgraphs == 1
    figure('units','normalized','innerposition',[0 0 1 1],'MenuBar','none')
    subplot(1,1,1) %1,3,1 (Line number,collumn number, graph id) - if you want to show more than 1 graph in same windows
    garbage.Xmax = 1500;
    garbage.Xmin = 0;
    garbage.Ymax = 1500;
    garbage.Ymin = 0;
    p = plot(Gaco,'XData',(ACOdataset.nodePosition(:,2)),'YData',(ACOdataset.nodePosition(:,3))); 
    line(ACOdataset.nodePosition(1:2,2),ACOdataset.nodePosition(1:2,3),'color','green','marker','o','linestyle','none','markersize',50)
    garbage.ax = gca;
    garbage.ax.XAxis.TickValues = 0:100:1000;
    garbage.ax.YAxis.TickValues = 0:100:1000;
    grid on
    hold on
    title(['Original WSN clone to ACO usage | ','Nodes number: ',num2str(ACOdataset.nodeNo),' | Nodes range: ', num2str(ACOdataset.range)])
    pause(2)
end

%% Main calcs section

%Get all neighbors from source an target
ACOneighborssource = neighbors(Gaco,ACOsourcenode);
ACOneighborstarget = neighbors(Gaco,ACOtargetnode);

%Get all source and target neighbors
counter=1;
for a = 1 : length(ACOneighborssource)
    for b = 1 : length(ACOneighborstarget)
        ACOallpossiblesourceANDtargets(counter,1)=ACOneighborssource(a,1);
        ACOallpossiblesourceANDtargets(counter,2)=ACOneighborstarget(b,1);
        counter=counter+1;
     end
end

%Get all Shortest paths between all source and target neighbors
for a = 1 : length(ACOallpossiblesourceANDtargets)
    ACOallpossibleshortestpaths(a).routes=shortestpathtree(Gaco,ACOallpossiblesourceANDtargets(a,1),ACOallpossiblesourceANDtargets(a,2));
end

%Get the total pheromone, energy and pathloss of each path
for a = 1 : length(ACOallpossibleshortestpaths)
    ACOtotalpathloss(a,:)=1/sum(ACOallpossibleshortestpaths(a).routes.Edges.Pathloss);
    ACOtotalenergy(a,:)=sum(ACOallpossibleshortestpaths(a).routes.Edges.Edgeenergy);
    ACOtotalpheromone(a,:)=ACOpheromone;
end

%First ACO call

for t = 1 : ACOmaxAntiteration
   
    %create ants
    colony=[];
    colony = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega);
    
end

% Fitness function - to know how good solutions found by colony are
for a = 1 : ACOantsNo
    tour = colony(a,:);
%     colony.ant(a).fitness = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour);
    fitness(a,1) = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone);
end

%Find the best fitness
ACOAllNeighborfitness = fitness(:,1);
[minVal,minIndex] = min(ACOAllNeighborfitness);
if minVal < ACObestfitness
    ACObestfitness = fitness(minIndex);
    ACObesttour=colony(minIndex);
end

    %Pheromone update
    ACOtotalpheromone(ACObesttour,1) = ACOtotalpheromone(ACObesttour,1)+ACOpheromone;

    %Evaporation event
    ACOtotalpheromone=(1 - ACOevapcoeficient).*ACOtotalpheromone;

    %Append first bebsttour index to best tour history
    ACObesttourhistory(end+1,1)=ACObesttour;
    
    
    %First call for nodes involved
    
    %Find all nodes involved in best tour event
    
    counter=1;
    for a = 1 : length(ACOallpossibleshortestpaths)
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
        ACOallnodesinvolvedinpossiblepaths(counter,1)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1);
        ACOallnodesinvolvedinpossiblepaths(counter,2)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
        counter=counter+1;
        end
    end
    ACOallnodesinvolvedinpossiblepaths=unique(ACOallnodesinvolvedinpossiblepaths);
        
    
    %Find nodes involved in besttour routing event
    ACObesttournodes=unique(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes); 
    
    %Import their best route nodes energy
    for a = 1 : length(ACObesttournodes)
       
        ACObesttournodes(a,2) = ACOdataset.nodePosition(ACObesttournodes(a,1),4);
        
    end
    
    %Find nodes NOT involved in besttour routing event
    ACOnodesNOTinvolvedinroutingevent=setdiff(ACOdataset.nodePosition(:,1),ACObesttournodes);
    
    
%     ACOnodesNOTinvolvedinroutingevent=setdiff(ACOallnodesinvolvedinpossiblepaths,ACObesttournodes);
%     for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
%            
%            ACOnodesNOTinvolvedinroutingevent(a,2)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4); 
%             
%     end
    
    for x = 1 : ACOenergycicles
    %Path energy usage
    %Involved nodes energy usage
    for a = 1 : length(ACObesttournodes)
        ACOdataset.nodePosition(ACObesttournodes(a,1),4)=ACOdataset.nodePosition(ACObesttournodes(a,1),4)-dataset.energyconsumptionperCicle^randomenergyfactor+dataset.energyrecoveryperCicle^randomenergyfactor;
    end
    
    for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
        if ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)<90
            ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)+dataset.energyrecoveryperCicle^randomenergyfactor;
            
        end
    end
    ACOpacket=ACOpacket+1;
            
    for a = 1 : length(ACOallpossibleshortestpaths)
        ACOtotalenergy(a,1)=0;
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
            ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
            ACOtotalenergy(a,1)=ACOtotalenergy(a,1)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1),4)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2),4);
        end
    end
    end

%Prepare output filename to save the reports
name=[strcat('WITH-ACO-simulation-',num2str(ACOdataset.nodeNo),'nodes-',datestr(now,'mm-dd-yyyy-HH-MM'),'.txt')];
fileID = fopen(fullfile(pwd,name),'w');
fprintf(fileID,'%6s %20s %20s %20s\r\n','|NodeNo|','|WITH ACO Scene|','|Packets sent|','|Best route|');



        
%% Main ACO loop and energy consumption

disp('Starting ants activity')

while ~isempty(ACOtotalenergy)
    
    ACOiterationcounter=ACOiterationcounter+1;
        
    if ACOdataset.nodePosition(:,4)<=0
        disp('Firt dead node detected')
        break
    end
    
    %Record best tour history
%     if ACObesttourhistory(end) ~= ACObesttour
        
    %Find all nodes involved in best tour event
    
    counter=1;
    for a = 1 : length(ACOallpossibleshortestpaths)
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
        ACOallnodesinvolvedinpossiblepaths(counter,1)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1);
        ACOallnodesinvolvedinpossiblepaths(counter,2)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
        counter=counter+1;
        end
    end
    ACOallnodesinvolvedinpossiblepaths=unique(ACOallnodesinvolvedinpossiblepaths);
        
    
    %Find nodes involved in besttour routing event
    ACObesttournodes=unique(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes); 
    
    %Import their best route nodes energy
    for a = 1 : length(ACObesttournodes)
       
        ACObesttournodes(a,2) = ACOdataset.nodePosition(ACObesttournodes(a,1),4);
        
    end
    
    %Find nodes NOT involved in besttour routing event
    ACOnodesNOTinvolvedinroutingevent=setdiff(ACOdataset.nodePosition(:,1),ACObesttournodes);
    
%     ACOnodesNOTinvolvedinroutingevent=setdiff(ACOallnodesinvolvedinpossiblepaths,ACObesttournodes);
    
    ACObesttourhistory(end+1,1)=ACObesttour;
    
%     end
for t = 1 : ACOenergycicles
    %Involved nodes energy usage
    for a = 1 : length(ACObesttournodes)
        ACOdataset.nodePosition(ACObesttournodes(a,1),4)=ACOdataset.nodePosition(ACObesttournodes(a,1),4)-dataset.energyconsumptionperCicle^randomenergyfactor+dataset.energyrecoveryperCicle^randomenergyfactor;
        
    end
    ACOpacket=ACOpacket+1;
    for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
        if ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)<90
            ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)+dataset.energyrecoveryperCicle^randomenergyfactor;
            
        end
    end
            
    for a = 1 : length(ACOallpossibleshortestpaths)
        ACOtotalenergy(a,1)=0;
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
            ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
            ACOtotalenergy(a,1)=ACOtotalenergy(a,1)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1),4)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2),4);
        end
    end
end


msg=['With ACO - Scene #: ',num2str(ACOiterationcounter),' | Packets sent: ', num2str(ACOpacket),' | Best route: ', num2str(ACObesttour),' | Routing nodes: ', num2str(ACObesttournodes(:,1)'),' | Nodes energy: ', num2str(ACObesttournodes(:,2)')];
%     msg=['With ACO - Scene #: ',num2str(ACOiterationcounter),' | Hops : ',num2str(length(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes)),' | Packets sent: ', num2str(ACOpacket),' | Dead node: ', num2str(deadnode),' | Routing nodes: ', num2str(ACObesttournodes)];
    disp(msg)

    
    report1=[dataset.nodeNo;iterationcounter;ACOpacket;ACObesttour];
    fprintf(fileID,'%6.0f %20.0f %20.0f %20.0f\r\n',report1);
    
    
    
    %Second ACO call

    for t = 1 : ACOmaxAntiteration

        %create ants
        colony=[];
        colony = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega);

    end

    % Fitness function - to know how good solutions found by colony are
    for a = 1 : ACOantsNo
        tour = colony(a,:);
    %     colony.ant(a).fitness = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour);
        fitness(a,1) = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone);
    end

    %Find the best fitness
    ACOAllNeighborfitness = fitness(:,1);
    [minVal,minIndex] = min(ACOAllNeighborfitness);
    if minVal < ACObestfitness
        ACObestfitness = fitness(minIndex);
        ACObesttour=colony(minIndex);
    end

        %Pheromone update
        ACOtotalpheromone(ACObesttour,1) = ACOtotalpheromone(ACObesttour,1)+ACOpheromone;

        %Evaporation event
        ACOtotalpheromone=(1 - ACOevapcoeficient).*ACOtotalpheromone;

    
    
    
end
        
fclose(fileID); %Close external data collector file


function [ colony ] = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega)

    for a = 1 : ACOantsNo
        for b = 1 : length(ACOtotalpheromone)
            
            ACOp_allNodes = ACOtotalpheromone(:,1).^ACOalpha.*ACOtotalenergy(:,1).^ACObeta; %.*ACOtotalpathloss(:,1).^ACOomega;

            P =  ACOp_allNodes./ sum(ACOp_allNodes);

            ACOnextroute = rouletteWheel(P);

            colony(a,b) = ACOnextroute;
%             colony2(b,:) = ACOnextroute;
%             colony.ant(a).tour = [colony.ant(a).tour, ACOnextroute];
         end

    end
end
%function [ fitness ] = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone)
% function [ fitness ] = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour)

%fitness = 0;

for i = 1 : length(tour)
    
    ACOcurrentroute = tour(i);
    %ACOnextnode = tour(i+1);
    
%     fitness = fitness + ACOedgesenergymatrix(1, ACOcurrentNode);
    fitness = fitness + ACOtotalpheromone(ACOcurrentroute, 1);% + (1/ACOtotalenergy(ACOcurrentroute, 1)); % + ACOtotalpathloss(ACOcurrentroute, 1);
%     fitness = fitness + ACOdesirebilitymatrix(ACOnextnode, 1);
    %fitness = fitness + Gaco.Edges(ACOcurrentNode, ACOnextnode);    

end


function [ ACOnextroute ] = rouletteWheel( P )

cumsumP = cumsum(P);

r = rand();

ACOnextroute = find(r <= cumsumP);

ACOnextroute = ACOnextroute(1);

end
%% ANOTA��ES IMPORTANTES