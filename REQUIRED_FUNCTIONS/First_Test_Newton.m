a=1:0.1:10;
ka_exp=1:0.01:10;
ka=10.^-ka_exp;
mu_max=1;
C_max=1;
C0=0.82;
x=1000;
y=0.87;
for i=1:size(a,2)
    for j=1:size(ka,2)
        mu_i(i,j)=mu_of_L_func(x,mu_max,Ka(j));
        y_i(i,j)=C_of_mu_func(mu_i(i,j),mu_max,C0,C_max,a(i));
        error_y(i,j)=(y-y_i(i,j))^2;
    end
end
[min_obj_fun,Indx_min_obj_fun]=min(error_y,[],'all','linear');
[Idx_a,Idx_ka]=ind2sub([size(error_y,1) size(error_y,2)],Indx_min_obj_fun);
min_a=a(Idx_a);
min_ka=ka(Idx_ka);
min_error=error_y(Idx_a,Idx_ka);
figure; hold on;
surf(error_y)
scatter3(Idx_ka,Idx_a,min_error,'filled','r');
