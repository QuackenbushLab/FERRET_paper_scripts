function A = run_grisli_single(varargin)
	file = varargin{1}
	start_spams;
	data_matrix = load(file);
	X = [transpose(data_matrix.ptime), transpose(data_matrix.X)];
	[~,I]=sort(X(:,1));
	X=X(I,:);
	Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);
	R=1500;
	L_array=20:10:90;
	alpha_min=.3;
	addpath(varargin{4});
	Rnk_array_TIGRESS_area_L=Compute_A_app_wo_ref(X,L_array,Alpha,R,alpha_min);
	A = Rnk_array_TIGRESS_area_L
	csvwrite(varargin{3}, A)
	return A;
end

