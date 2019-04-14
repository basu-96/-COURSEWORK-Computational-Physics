n_sites = 2;
e1 = 0.0;
e2 = 4.0;
t = 1.0;
ham_mat(1,1) = e1;
ham_mat(1,2) = -t;
ham_mat(2,1) = -t;
ham_mat(2,2) = e2;
eta = 0.01;
omega = [-2*t:0.01:2*t];
plotdata = [];
for i = 1:numel(omega)
	z = omega(i) + 1i*eta;
	g11 = 1/(z - e1);
	g22 = 1/(z - e2);
	G11 = 1/((1/g11 - t*g22*t));
	G22 = 1/((1/g22 - t*g11*t));
	-imag(G11)
	% plotdata = [plotdata;[omega -imag(G11)/pi -imag(G22)/pi]];
end
% figure;
% plot(plotdata(:,1), plotdata(:,2));
% saveas(gcf, 'G11.png');
