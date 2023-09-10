Nsites=10;%Lattice scale for one dimensional lattice

Str_single={};
%Generate two site operator H=sum_{i}J1*b_{i}^{diagger}b_{i+1}^{minus} (nearest hopping)
hopping1={'J1*b_diagger','b_minus',1};
for i=1:Nsites-hopping1{3}
    j=i+hopping1{3};
    first_str=hopping1{1};second_str=hopping1{2};
    row_end=size(Str_single,1);
    for loop=1:Nsites
        if loop==i
            Str_single{row_end+1,loop}=first_str;
        elseif loop==j
            Str_single{row_end+1,loop}=second_str;
        else
            Str_single{row_end+1,loop}='I';
        end
    end

    j=i+hopping1{3};
    first_str=hopping1{2};second_str=hopping1{1};
    row_end=size(Str_single,1);
    for loop=1:Nsites
        if loop==i
            Str_single{row_end+1,loop}=first_str;
        elseif loop==j
            Str_single{row_end+1,loop}=second_str;
        else
            Str_single{row_end+1,loop}='I';
        end
    end
end

%Generate two site operator H=sum_{i}J2*n_{i}n_{i+1} (long range interaction)
% hopping1={'U2*n','n',1};
% for i=1:Nsites-hopping1{3}
%     j=i+hopping1{3};
%     first_str=hopping1{1};second_str=hopping1{2};
%     row_end=size(Str_single,1);
%     for loop=1:Nsites
%         if loop==i
%             Str_single{row_end+1,loop}=first_str;
%         elseif loop==j
%             Str_single{row_end+1,loop}=second_str;
%         else
%             Str_single{row_end+1,loop}='I';
%         end
%     end
% 
%     j=i+hopping1{3};
%     first_str=hopping1{2};second_str=hopping1{1};
%     row_end=size(Str_single,1);
%     for loop=1:Nsites
%         if loop==i
%             Str_single{row_end+1,loop}=first_str;
%         elseif loop==j
%             Str_single{row_end+1,loop}=second_str;
%         else
%             Str_single{row_end+1,loop}='I';
%         end
%     end
% end

%Generate two site operator H=sum_{i}J2*b_{i}^{diagger}b_{i+2}^{minus} (long range hopping)
hopping1={'J2*b_diagger','b_minus',2};
for i=1:Nsites-hopping1{3}
    j=i+hopping1{3};
    first_str=hopping1{1};second_str=hopping1{2};
    row_end=size(Str_single,1);
    for loop=1:Nsites
        if loop==i
            Str_single{row_end+1,loop}=first_str;
        elseif loop==j
            Str_single{row_end+1,loop}=second_str;
        else
            Str_single{row_end+1,loop}='I';
        end
    end

    j=i+hopping1{3};
    first_str=hopping1{2};second_str=hopping1{1};
    row_end=size(Str_single,1);
    for loop=1:Nsites
        if loop==i
            Str_single{row_end+1,loop}=first_str;
        elseif loop==j
            Str_single{row_end+1,loop}=second_str;
        else
            Str_single{row_end+1,loop}='I';
        end
    end
end

%Generate three site operator H=sum_{i}Omega*P_{i}X_{i+1}P_{i+2} (PXP model)
% hopping2={'Omega*P','X','P',1,2};
% for i=1:Nsites-hopping2{5}
%     j=i+hopping2{4};j2=i+hopping2{5};
%     first_str=hopping2{1};second_str=hopping2{2};Third_str=hopping2{3};
%     row_end=size(Str_single,1);
%     for loop=1:Nsites
%         if loop==i
%             Str_single{row_end+1,loop}=first_str;
%         elseif loop==j
%             Str_single{row_end+1,loop}=second_str;
%         elseif loop==j2
%             Str_single{row_end+1,loop}=Third_str;
%         else
%             Str_single{row_end+1,loop}='I';
%         end
%     end
% end

%Generate onsite operator H0=mu*n+U*n(n-1)/2(on site potention plus interaction)
onsite={'H0'};
for i=1:Nsites
    row_end=size(Str_single,1);
    for loop=1:Nsites
        if loop==i
            Str_single{row_end+1,loop}=onsite{1};
        else
            Str_single{row_end+1,loop}='I';
        end
    end
end
%%
%Main program
left_number=1;left_index={1:size(Str_single,1)};
Total_MPO={};
for loop=1:Nsites-1
    new_left_number=0;
    new_left_index={};
    Mid_A={};
    for kk=1:left_number
        A=categorical(Str_single(left_index{kk},loop));
        category=categories(A);
        for mm=1:length(category)
            new_left_number=new_left_number+1;
            Mid_A{kk,new_left_number}=category{mm};
            site=strcmp(Str_single(left_index{kk},loop),category{mm});
            new_left_index{new_left_number}=left_index{kk}(site==1);
        end
    end
    Total_MPO{loop}=Mid_A;
    left_number=new_left_number;
    left_index=new_left_index;
end
Mid_A={};
for kk=1:left_number
    A=categorical(Str_single(left_index{kk},Nsites));
    category=categories(A);
    for mm=1:length(category)
        Mid_A{kk,1}=category{mm};
    end
end
Total_MPO{Nsites}=Mid_A;

right_number=1;right_index={1};
for loop=Nsites:-1:2
    Mid_A=Total_MPO{loop};
    new_right_number=0;
    new_right_index={};
    new_A={};
    for kk=1:right_number
        Str_single_site=Mid_A(:,right_index{kk});
        [index_matrix] = Find_empyty(Str_single_site);
        A=categorical(Str_single_site(index_matrix==1));
        category=categories(A);
        for mm=1:length(category)
            site=strcmp(Str_single_site,category{mm});
            [x,y]=find(site==1);
            new_A{min(x),kk}=category{mm};
            new_right_index{min(x)}=x;
        end
    end
    [index_matrix] = Find_empyty(new_A);
    row_index=sum(index_matrix,2);
    Total_MPO{loop}=new_A(row_index>0,:);
    right_number=size(Total_MPO{loop},1);
    right_index=new_right_index(row_index>0);
end

%%

for loop=1:Nsites
    Total_MPO{loop}
end