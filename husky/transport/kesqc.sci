function [TE] = computeTE(Hs,Hc,Hi,Hi1,Hi2,energies)

		p=1;
		
		for e=energies
		
		// on forme le propagateur
		// Attention pas de classement => electrodes de dim 1
		P=[
		zeros(Hc) , eye(Hc)
		-inv(Hi)*Hi' , -inv(Hi)*(Hc-e*eye(Hc))];
		
		// on trouve les vecteurs propre du propagateur
		[x,h]=spec(P);
		
		// on forme la matrice du systeme
		S=[
		Hs-e*eye(Hs)  Hi1          zeros(Hi1)   Hi2          zeros(Hi2)   zeros(Hi1) zeros(Hi2)
		Hi1'          Hc-e*eye(Hc) Hi           zeros(Hc)    zeros(Hc)    zeros(Hc)  zeros(Hc)
		zeros(Hi1')   Hi'          Hc-e*eye(Hc) zeros(Hc)    zeros(Hc)    Hi         zeros(Hc)
		Hi2'          zeros(Hc)    zeros(Hc)    Hc-e*eye(Hc) Hi           zeros(Hc)  zeros(Hc)
		zeros(Hi2')   zeros(Hc)    zeros(Hc)    Hi'          Hc-e*eye(Hc) zeros(Hc)  Hi
		];
		
		// on trouve le noyau de S
		Ker=kernel(S);
		
		// on extrait les composantes sur les cellules 2 et 3 des electrodes
		V=[
		Ker(size(Hs,1)+size(Hc,1)+[1:size(Hc,1)],:)
		Ker(size(Hs,1)+4*size(Hc,1)+[1:size(Hc,1)],:)
		Ker(size(Hs,1)+3*size(Hc,1)+[1:size(Hc,1)],:)
		Ker(size(Hs,1)+5*size(Hc,1)+[1:size(Hc,1)],:)
		];
		
		// passe dans la base des etats entrants et sortants
		X=[
		x(:,1)        zeros(x(:,1))  x(:,2)        zeros(x(:,2))
		zeros(x(:,1)) x(:,1)         zeros(x(:,2)) x(:,2)
		];
		V=inv(X)*V;
		
		// on recupere les coefficients des entrants et sortants
		M1=V([1:2*size(Hc,1)],:);
		M2=V(2*size(Hc,1)+[1:2*size(Hc,1)],:);
		
		// on calcule la matrice de Scatering
// 		if size(M1,1)~=size(M1,2)
// 			printf("pro
// 			pause;
// 			
// 		end
		Scat=M2*inv(M1);
	
		// Attention pas de renormalisation ici
		// on passe directement dans la base des courants
		T=abs(Scat).^2;
		
		
		TE(p)=T(1,2);
			
		p=p+1;

  end

endfunction


function [TE] = calculTE(H,energies,n1,n2,eps,e0,h)

	//////////////////////////////////////
	//	definiton des electrodes			
	//////////////////////////////////////
	
	// hamiltonien de la cellule de l'electrode
	Hc=[
	e0	
	];
	
	// hamiltonien d'interaction des cellules
	// de l'electrode (n+1) => n
	Hi=[
	h	
	];
			
	////////////////////////////////////////
	//	interaction		
	////////////////////////////////////// 
	
	//interaction electrode 1 - systeme
	Hi1 = zeros(1,size(H,1));
	Hi1(1,n1) = eps;

	
	// Hi1 = U'*Hi1;
	
	// interaction electrode 2 - systeme
	Hi2 = zeros(1,size(H,1));
	Hi2(1,n2) = eps;

	[u,V]=spec(H);
	Hi1_=u*Hi1';
	Hi2_=u*Hi2';


	//////////////////////////////////////
	//	calcul du te a partir des TE			
	//////////////////////////////////////

	TE = computeTE(H,Hc,Hi,Hi1',Hi2',energies);


endfunction



function [TE] = calculTE_assym(H,energies,eps1,eps2,e0,h)

	//////////////////////////////////////
	//	definiton des electrodes			
	//////////////////////////////////////
	
	// hamiltonien de la cellule de l'electrode
	Hc=[
	e0	
	];
	
	// hamiltonien d'interaction des cellules
	// de l'electrode (n+1) => n
	Hi=[
	h	
	];
			
	////////////////////////////////////////
	//	interaction		
	////////////////////////////////////// 
	
	//interaction electrode 1 - systeme
	Hi1 = zeros(1,size(H,1));
	Hi1 = eps1;

	
	// Hi1 = U'*Hi1;
	
	// interaction electrode 2 - systeme
	Hi2 = zeros(1,size(H,1));
	Hi2 = eps2;
	
	//////////////////////////////////////
	//	calcul du te a partir des TE			
	//////////////////////////////////////

	TE = computeTE(H,Hc,Hi,Hi1',Hi2',energies);


endfunction



function [TE] = calculTE_3elec(H,energies,n1,n2,n3,eps,e0,h)

	//////////////////////////////////////
	//	definiton des electrodes			
	//////////////////////////////////////
	
	// hamiltonien de la cellule de l'electrode
	Hc=[
	e0	
	];
	
	// hamiltonien d'interaction des cellules
	// de l'electrode (n+1) => n
	Hi=[
	h	
	];
			
	////////////////////////////////////////
	//	interaction		
	////////////////////////////////////// 
	
	//interaction electrode 1 - systeme
	Hi1 = zeros(1,size(H,1));
	Hi1(1,n1) = eps;

	
	// interaction electrode 2 - systeme
	Hi2 = zeros(1,size(H,1));
	Hi2(1,n2) = eps;


	// interaction electrode 2 - systeme
	Hi3 = zeros(1,size(H,1));
	Hi3(1,n3) = eps;


	//////////////////////////////////////
	//	calcul du te a partir des TE			
	//////////////////////////////////////
	
	TE=computeTE_3el(H,Hi1',Hi2',Hi3',Hc,Hi,energies);


endfunction

function [TE] = calculTE_4elec(H,energies,n1,n2,n3,n4,eps,e0,h)

	//////////////////////////////////////
	//	definiton des electrodes			
	//////////////////////////////////////
	
	// hamiltonien de la cellule de l'electrode
	Hc=[
	e0	
	];
	
	// hamiltonien d'interaction des cellules
	// de l'electrode (n+1) => n
	Hi=[
	h	
	];
			
	////////////////////////////////////////
	//	interaction		
	////////////////////////////////////// 
	
	//interaction electrode 1 - systeme
	Hi1 = zeros(1,size(H,1));
	Hi1(1,n1) = eps;

	
	// interaction electrode 2 - systeme
	Hi2 = zeros(1,size(H,1));
	Hi2(1,n2) = eps;


	// interaction electrode 3 - systeme
	Hi3 = zeros(1,size(H,1));
	Hi3(1,n3) = eps;
	
	// interaction electrode 4 - systeme
	Hi4 = zeros(1,size(H,1));
	Hi4(1,n4) = eps;



	//////////////////////////////////////
	//	calcul du te a partir des TE			
	//////////////////////////////////////
	
	TE=computeTE_4el(H,Hi1',Hi2',Hi3',Hi4',Hc,Hi,energies);


endfunction


