function [out]=loadFEM(ii)

    if ii>9
        error('estendere funzione per ii>9')
    end

    num_file = num2str(ii);

    str = strcat('diagn_','aa','00',num_file,'.csv');
    aa = load(str);
    aa = aa';
    out.aa = complex(aa(:,1),aa(:,2));
    
    str = strcat('diagn_','tn','00',num_file,'.csv');
    tn = load(str);
    tn = tn';
    out.tn = complex(tn(:,1),tn(:,2));
    
    str = strcat('diagn_','sol','00',num_file,'.csv');
    sol = load(str); 
    sol = sol';
    out.sol = complex(sol(:,1),sol(:,2));
    
    str = strcat('diagn_','s','00',num_file,'.csv');
    s = load(str);
    out.s = complex(s(1),s(2));
    
    str = strcat('diagn_','current','00',num_file,'.csv');
    I = load(str);
    out.I = complex(I(:,1),I(:,2));
end