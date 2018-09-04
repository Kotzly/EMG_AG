
offset=0;

% Saves the movements from the NinaPro .mat files to individual .csv files
for ex = 1:3
    
    exName=strcat('E',num2str(ex));
    load(strcat('S1_',exName,'_A1.mat'))
    
    for i = 11:12
        if i == 11
            muscle='bic';
        else
            muscle='tri';
        end
        temp=[];
        for j= 2:size(emg,1)

            if restimulus(j)~=0


                if restimulus(j)+offset<=9
                    o='0';
                else
                    o='';
                end
                movement=strcat(o,num2str(restimulus(j)+offset));
                rep=num2str(rerepetition(j));
                
%                arqName=strcat(muscle,'_',exName,'_',movement,'_',rep,'.csv');
                arqName=strcat(muscle,'_',movement,'_',rep,'.csv');
                temp(end+1)=emg(j,i);
            elseif restimulus(j-1)~=0
                if ~ exist(arqName,'file')
                    arqName
                    csvwrite(strcat('csv\',arqName),temp);
                end
                temp=[];
            end
        end
    end
    
    offset=offset+max(stimulus);
end