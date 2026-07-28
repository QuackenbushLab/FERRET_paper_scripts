function A = run_grisli_single(varargin)
	file = varargin{1}
	fprintf(file)
	start_spams;
	data_matrix = load(file);
	fprintf(data_matrix)
	X = [transpose(data_matrix.ptime), transpose(data_matrix.X)];
	[~,I]=sort(X(:,1));
	X=X(I,:);
	Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);
	R=1500;
	L_array=20:10:90;
	alpha_min=.3;
	addpath('/home/ubuntu/GRISLI/GRISLI/');
	Rnk_array_TIGRESS_area_L=Compute_A_app_wo_ref(X,L_array,Alpha,R,alpha_min);
	A = Rnk_array_TIGRESS_area_L
	csvwrite('/home/ubuntu/GRISLI_HTAN/file', A)
	return A;
end

