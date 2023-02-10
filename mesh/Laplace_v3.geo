SetFactory("OpenCASCADE");

Macro Def_physical_curve
	Physical Curve(Str(name_physical_curve),i_physical_curve) = raccolta_curve();
Return

Macro Def_physical_surface
	Physical Surface(Str(name_physical_surface),i_physical_surface) = raccolta_superfici();
Return

Macro GenCirControllo // circonferenza di controllo
	circControllo = newl;
	Circle(circControllo) = {x, y, 0, R_circControllo, 0, 2*Pi};

	If(y > 0)
		Curve{circControllo} In Surface{S_pl};
	ElseIf(y < 0)
		Curve{circControllo} In Surface{S_coppa_p};
	EndIf
	
	Transfinite Curve {circControllo} = n_nodi_cirControllo Using Progression 1;
Return

Macro GenEllControllo // ellisse di controllo
	circControllo = newl;
	// Printf("newl = %g",circControllo);
	Ellipse(circControllo) = {x, y, 0, 0.55*L_x, 10*L_y, 0, 2*Pi};

	If(y > 0)
		Curve{circControllo} In Surface{S_pl};
	ElseIf(y < 0)
		Curve{circControllo} In Surface{S_coppa_p};
	EndIf
	
	Transfinite Curve {circControllo} = n_nodi_cirControllo Using Progression 1;
Return

Macro GenCondCil // conduttore cilindrico
	AA = newp;
	Point(AA) = {x, y, 0, lms};
	A = newp;
	Point(A) = {x,y+R_ext, 0, lms};
	B = newp;
	Point(B) = {x-R_ext, y, 0, lms};
	C = newp;
	Point(C) = {x, y-R_ext, 0, lms};
	D = newp;
	Point(D) = {x+R_ext, y, 0, lms};
	E = newp;
	Point(E) = {x, y+R_int, 0, lms};
	F = newp;
	Point(F) = {x-R_int, y, 0, lms};
	G = newp;
	Point(G) = {x, y-R_int, 0, lms};
	H = newp;
	Point(H) = {x+R_int, y, 0, lms};

	// contorno esterno cond di fase
	Cir_A =  newl;
	Circle(Cir_A) = {A, AA, B};
	Cir_B =  newl;
	Circle(Cir_B) = {B, AA, C};
	Cir_C =  newl;
	Circle(Cir_C) = {C, AA, D};
	Cir_D =  newl;
	Circle(Cir_D) = {D, AA, A};

	// contorno interno cond di fase
	Cir_E =  newl;
	Circle(Cir_E) = {E, AA, F};
	Cir_F =  newl;
	Circle(Cir_F) = {F, AA, G};
	Cir_G =  newl;
	Circle(Cir_G) = {G, AA, H};
	Cir_H =  newl;
	Circle(Cir_H) = {H, AA, E};
	l_AE = newl;
	Line(l_AE) = {A, E};
	l_BF = newl;
	Line(l_BF) = {B, F};
	l_CG = newl;
	Line(l_CG) = {C, G};
	l_DH = newl;
	Line(l_DH) = {D, H};

	// 4 sup esterne condutt fase, saranno strutt
	S_1 = news;
	ll_S1 = newll;
	Curve Loop(ll_S1) = {Cir_A, l_BF, -Cir_E, -l_AE};
	Plane Surface(S_1) = {ll_S1};
	S_2 = news;
	ll_S2 = newll;
	Curve Loop(ll_S2) = {Cir_B, l_CG, -Cir_F, -l_BF};
	Plane Surface(S_2) = {ll_S2};
	S_3 = news;
	ll_S3 = newll;
	Curve Loop(ll_S3) = {Cir_C, l_DH, -Cir_G, -l_CG};
	Plane Surface(S_3) = {ll_S3};
	S_4 = news;
	ll_S4 = newll;
	Curve Loop(ll_S4) = {Cir_D, l_AE, -Cir_H, -l_DH};
	Plane Surface(S_4) = {ll_S4};
	// sup centrale, non strutturata
	S_centro = news;
	ll_S_centro = newll;
	Curve Loop(ll_S_centro) = {Cir_E, Cir_F, Cir_G, Cir_H};
	Plane Surface(S_centro) = {ll_S_centro};
	Point{AA} In Surface{S_centro};

	Transfinite Surface {S_1};
	Transfinite Surface {S_2};
	Transfinite Surface {S_3};
	Transfinite Surface {S_4};

	Transfinite Curve {Cir_A, Cir_B, Cir_C, Cir_D} = n_nodi_azimut Using Progression 1;
	Transfinite Curve {Cir_E, Cir_F, Cir_G, Cir_H} = n_nodi_azimut Using Progression 1;

	Transfinite Curve {l_AE, l_BF, l_CG, l_DH} = n_nodi_radial Using Progression progr_radial;

	If(y > 0)
		bb() = BooleanDifference{Surface{S_pl}; Delete;}{Surface{S_1, S_2, S_3, S_4, S_centro};};
	ElseIf(y < 0)
		bb() = BooleanDifference{Surface{S_coppa_p}; Delete;}{Surface{S_1, S_2, S_3, S_4, S_centro};};
	EndIf
	
	If(flag_isPipe == 0)
		raccolta_superfici = {S_1, S_2, S_3, S_4, S_centro};
		Call Def_physical_surface;
	ElseIf(flag_isPipe == 1)
		name_physical_surface = "pipe";
		raccolta_superfici = {S_1, S_2, S_3, S_4}; // PIPELINE
		Call Def_physical_surface;
		
		name_physical_surface = "pipe_centro";
		raccolta_superfici = {S_centro}; // CENTRO PIPELINE
		i_physical_surface = i_physical_surface+1;
		Call Def_physical_surface;
	EndIf

Return

Macro GenCondRectUnstruct // conduttore rettangolare

	lms = 5E-4;
	
	p_cntr = newp;
	Point(p_cntr) = {x, y, 0, lms*2};

	p_A = newp;
	Point(p_A) = {x-L_x/2,y-L_y/2, 0, lms};
	p_B = newp;
	Point(p_B) = {x-L_x/2, y-L_y/2+L_y, 0, lms};
	p_C = newp;
	Point(p_C) = {x-L_x/2+L_x, y-L_y/2+L_y, 0, lms};
	p_D = newp;
	Point(p_D) = {x-L_x/2+L_x, y-L_y/2, 0, lms};

	If (flag_rotate_mit == 1)
        Rotate{{0,0,1}, {x, y,0}, Pi/4}{ Point{p_A,p_B,p_C,p_D}; } // ruoto attorno ad asse Z, using come centro di rotazione il punto di coord (x,y)
    EndIf
	
	// linee contorno ext
	l_AB = newl;
	Line(l_AB) = {p_A, p_B};
	l_BC = newl;
	Line(l_BC) = {p_B, p_C};
	l_CD = newl;
	Line(l_CD) = {p_C, p_D};
	l_DA = newl;
	Line(l_DA) = {p_D, p_A};

	// linee 4 spicchi interni
	l_Acntr = newl;
	Line(l_Acntr) = {p_A, p_cntr};
	l_Bcntr = newl;
	Line(l_Bcntr) = {p_B, p_cntr};
	l_Ccntr = newl;
	Line(l_Ccntr) = {p_C, p_cntr};
	l_Dcntr = newl;
	Line(l_Dcntr) = {p_D, p_cntr};
	
	/*
	ll_S1 = newll;
	Curve Loop(ll_S1) = {l_AB, l_BC, l_CD, l_DA};
	S_1 = news;
	Plane Surface(S_1) = {ll_S1};
	If(y > 0)
		BooleanDifference{ Surface{S_air}; Delete; }{ Surface{S_1};}
	ElseIf(y < 0)
		BooleanDifference{ Surface{S_soil}; Delete; }{ Surface{S_1};}
	EndIf

	Physical Surface(i_physical_surface) = {S_1};
	*/
	
	/*
	Transfinite Curve {l_AB} = n_nodi_radial Using Progression progr_radial;
	Transfinite Curve {l_BC} = n_nodi_radial Using Progression progr_radial;
	Transfinite Curve {l_CD} = n_nodi_radial Using Progression progr_radial;
	Transfinite Curve {l_DA} = n_nodi_radial Using Progression progr_radial;
	*/
	
	ll_S1 = newll;
	Curve Loop(ll_S1) = {l_Acntr, l_AB, l_Bcntr};
	S_1 = news;
	Plane Surface(S_1) = {ll_S1};
	// Transfinite Surface {S_1};

	ll_S2 = newll;
	Curve Loop(ll_S2) = {l_Bcntr, l_BC, l_Ccntr};
	S_2 = news;
	Plane Surface(S_2) = {ll_S2};
	// Transfinite Surface {S_2};

	ll_S3 = newll;
	Curve Loop(ll_S3) = {l_Ccntr, l_CD, l_Dcntr};
	S_3 = news;
	Plane Surface(S_3) = {ll_S3};
	// Transfinite Surface {S_3};

	ll_S4 = newll;
	Curve Loop(ll_S4) = {l_Dcntr, l_DA, l_Acntr};
	S_4 = news;
	Plane Surface(S_4) = {ll_S4};
	// Transfinite Surface {S_4};

	If(y > 0)
		BooleanDifference{ Surface{S_pl}; Delete; }{ Surface{S_1,S_2,S_3,S_4};}
	ElseIf(y < 0)
		BooleanDifference{ Surface{S_coppa_p}; Delete; }{ Surface{S_1,S_2,S_3,S_4};}
	EndIf

	Physical Surface(i_physical_surface) = {S_1,S_2,S_3,S_4};
	
	/*
	If(x > 0)
		pto_segue = newp;
		Point(pto_segue) = {x,-1E-3, 0, 0.01};
		Point{pto_segue} In Surface{S_soil};
	EndIf
	*/

Return

Macro GenCondRectStruct // conduttore rettangolare

	lms = 1;

	p_A = newp;
	Point(p_A) = {x-L_x/2,y-L_y/2, 0, lms};
	p_B = newp;
	Point(p_B) = {x-L_x/2, y-L_y/2+L_y, 0, lms};
	p_C = newp;
	Point(p_C) = {x-L_x/2+L_x, y-L_y/2+L_y, 0, lms};
	p_D = newp;
	Point(p_D) = {x-L_x/2+L_x, y-L_y/2, 0, lms};

	If (flag_rotate_mit == 1)
        Rotate{{0,0,1}, {x, y,0}, Pi/4}{ Point{p_A,p_B,p_C,p_D}; } // ruoto attorno ad asse Z, using come centro di rotazione il punto di coord (x,y)
    EndIf
	
	// linee contorno ext
	l_AB = newl;
	Line(l_AB) = {p_A, p_B};
	l_BC = newl;
	Line(l_BC) = {p_B, p_C};
	l_CD = newl;
	Line(l_CD) = {p_C, p_D};
	l_DA = newl;
	Line(l_DA) = {p_D, p_A};
	
	aspect_ratio = L_x/L_y;
	
	aaa = 1.05;
	bbb = 1.0;
	n_p_s_x = 20;
	n_p_s_y = 15;
	
	Transfinite Curve {l_AB} = n_p_s_y Using Progression 1/aaa;
	Transfinite Curve {l_BC} = n_p_s_x Using Progression bbb;
	Transfinite Curve {l_CD} = n_p_s_y Using Progression aaa;
	Transfinite Curve {l_DA} = n_p_s_x Using Progression bbb;
	
	ll_S1 = newll;
	Curve Loop(ll_S1) = {l_AB, l_BC, l_CD, l_DA};
	S_1 = news;
	Plane Surface(S_1) = {ll_S1};
	If(y > 0)
		BooleanDifference{ Surface{S_pl}; Delete; }{ Surface{S_1};}
	ElseIf(y < 0)
		BooleanDifference{ Surface{S_coppa_p}; Delete; }{ Surface{S_1};}
	EndIf
	
	Transfinite Surface {S_1};
	Physical Surface(i_physical_surface) = {S_1};
	
Return

Macro GenStripRect // lastra
	lms = 1;

	p_A = newp;
	Point(p_A) = {x-L_x/2,y-L_y/2, 0, lms};
	p_B = newp;
	Point(p_B) = {x-L_x/2, y-L_y/2+L_y, 0, lms};
	p_C = newp;
	Point(p_C) = {x-L_x/2+L_x, y-L_y/2+L_y, 0, lms};
	p_D = newp;
	Point(p_D) = {x-L_x/2+L_x, y-L_y/2, 0, lms};

	If (flag_rotate_mit == 1)
        Rotate{{0,0,1}, {x, y,0}, Pi/4}{ Point{p_A,p_B,p_C,p_D}; } // ruoto attorno ad asse Z, using come centro di rotazione il punto di coord (x,y)
    EndIf
	
	// linee contorno ext
	l_AB = newl;
	Line(l_AB) = {p_A, p_B};
	l_BC = newl;
	Line(l_BC) = {p_B, p_C};
	l_CD = newl;
	Line(l_CD) = {p_C, p_D};
	l_DA = newl;
	Line(l_DA) = {p_D, p_A};
	
	aspect_ratio = L_x/L_y;
	
	aaa = 1.05;
	bbb = 1.0;
	n_p_s_x = 250;
	n_p_s_y = 15;
	
	Transfinite Curve {l_AB} = n_p_s_y Using Progression 1/aaa;
	Transfinite Curve {l_BC} = n_p_s_x Using Progression bbb;
	Transfinite Curve {l_CD} = n_p_s_y Using Progression aaa;
	Transfinite Curve {l_DA} = n_p_s_x Using Progression bbb;
	
	ll_S1 = newll;
	Curve Loop(ll_S1) = {l_AB, l_BC, l_CD, l_DA};
	S_1 = news;
	Plane Surface(S_1) = {ll_S1};
	If(y > 0)
		BooleanDifference{ Surface{S_pl}; Delete; }{ Surface{S_1};}
	ElseIf(y < 0)
		BooleanDifference{ Surface{S_coppa_p}; Delete; }{ Surface{S_1};}
	EndIf
	
	Transfinite Surface {S_1};
	Physical Surface(i_physical_surface) = {S_1};
	
Return

// ******************************************************************* //
// ******************************************************************* //

// IMPOSTAZIONI PIPE
n_spost = 1; // 15
x_pipe_racc(0) = 2;
/*
x_pipe_racc(1) = 3;
x_pipe_racc(2) = 0;
x_pipe_racc(3) = 3;
x_pipe_racc(4) = 0;
x_pipe_racc(5) = 3;
x_pipe_racc(6) = 6;
x_pipe_racc(7) = 10;
x_pipe_racc(8) = 15;
x_pipe_racc(9) = 20;
x_pipe_racc(10) = 25;
x_pipe_racc(11) = 30;
x_pipe_racc(12) = 35;
x_pipe_racc(13) = 40;
x_pipe_racc(14) = 45;
*/

y_pipe = -1.1;
R_ext_pipe = 0.5;		
R_int_pipe = R_ext_pipe - 0.015;

// IMPOSTAZIONI COND FASE
// spost orizz terne
x_t_1 = 0;
h_base_traliccio = 0; // 28.45 + 13.25
R_ext_fase = 0.015;
R_int_fase = R_ext_fase/2;

xx[0] = 3.5 + x_t_1; // fase 1
yy[0] = h_base_traliccio + 12;

xx[1] = -3 + x_t_1; // fase sx
yy[1] = h_base_traliccio + 14;

xx[2] = 2.9 + x_t_1; // fase dx
yy[2] = h_base_traliccio + 16;

// IMPOSTAZIONI OGW
xx[3] = 0 + x_t_1; // OGW SX
yy[3] = h_base_traliccio + 20;

R_ext_OGW = 0.006;
R_int_OGW = R_ext_OGW/2;

// IMPOSTAZIONI TERRENO
x_0 = 0;
y_0 = 0;
R_ext_domain = 6E2; // 10E3
dim_pti_dominio_ext = 50.0; // 

// IMPOSTAZIONI MIT
altezza_mit = -0.25;

R_mit = 0.008;
lms_mit = 0.0001;

// MAIN LOOP
For i_spost In {0:n_spost-1}

	x_pipe = x_pipe_racc(i_spost);
	
	n_mit = 0;
	xx_mit[n_mit] = x_pipe-1.6;	// mit 1, sx
	yy_mit[n_mit] = altezza_mit;
	n_mit ++;

	xx_mit[n_mit] = x_pipe-0.55; // mit 2
	yy_mit[n_mit] = altezza_mit;
	n_mit ++;

	xx_mit[n_mit] = x_pipe+0.55; // mit 3
	yy_mit[n_mit] = altezza_mit;
	n_mit ++;

	xx_mit[n_mit] = x_pipe+1.6;	// mit 4
	yy_mit[n_mit] = altezza_mit;
	n_mit ++;

	// creo i 4 punti esterni del dominio
	A_ext = newp;
	Point(newp) = {R_ext_domain, 0, 0, dim_pti_dominio_ext};
	B_ext = newp;
	Point(newp) = {0, R_ext_domain, 0, dim_pti_dominio_ext};
	C_ext = newp;
	Point(newp) = {-R_ext_domain, 0, 0, dim_pti_dominio_ext};
	D_ext = newp;
	Point(newp) = {0, -R_ext_domain, 0, dim_pti_dominio_ext};

	// punto centrale
	centro_ext = newp;
	Point(centro_ext) = {x_0+0, y_0+0, 0, 0.5};

	// punto al livello del terreno che si sposta con pipeline per garantire dimensione mesh
	pto_segue_pipe = newp;
	Point(pto_segue_pipe) = {x_pipe+0, y_0+0, 0, 0.01};

	// creo archi di cerchio per connettere punti esterni del dominio
	Cir_A = newl;
	Circle(Cir_A) = {A_ext, centro_ext, B_ext};
	Cir_B = newl;
	Circle(Cir_B) = {B_ext, centro_ext, C_ext};
	Cir_C = newl;
	Circle(Cir_C) = {C_ext, centro_ext, D_ext};
	Cir_D = newl;
	Circle(Cir_D) = {D_ext, centro_ext, A_ext};

    // cancello punto centro geometria, solo per costruzione
    Recursive Delete {
        Point{centro_ext}; 
    }

    // coppa che contiene mitigation e pipe
    R_coppa_pipe = 10; // R_ext_pipe*5
    lms_coppa_pipe = 0.1;
	sotto = newp; Point(sotto) = {x_pipe, -R_coppa_pipe+y_pipe, 0, lms_coppa_pipe};
    p_coppa_p_dx = newp; Point(p_coppa_p_dx) = {x_pipe+R_coppa_pipe-y_pipe, 0, 0, lms_coppa_pipe};
	p_coppa_p_sx = newp; Point(p_coppa_p_sx) = {x_pipe-R_coppa_pipe+y_pipe, 0, 0, lms_coppa_pipe};
    semic_1 = newl; Circle(semic_1) = {p_coppa_p_sx, pto_segue_pipe, sotto};
	semic_2 = newl; Circle(semic_2) = {sotto, pto_segue_pipe, p_coppa_p_dx};

    l_coppa_p_sx = newl; Line(l_coppa_p_sx) = {pto_segue_pipe,p_coppa_p_sx};
    l_coppa_p_dx = newl; Line(l_coppa_p_dx) = {pto_segue_pipe,p_coppa_p_dx};

    Transfinite Curve {l_coppa_p_sx} = 200 Using Progression 1.02;
    Transfinite Curve {l_coppa_p_dx} = 200 Using Progression 1.02;

    Transfinite Curve {semic_1} = 50 Using Progression 1.0;
    Transfinite Curve {semic_2} = 50 Using Progression 1.0;

    ll_int_S_air[0] = newll;
    Curve Loop(ll_int_S_air[0]) = {l_coppa_p_sx, semic_1, semic_2, l_coppa_p_dx};

    // SUPERFICIE COPPA PIPELINE
    S_coppa_p = news; Plane Surface(S_coppa_p) = {ll_int_S_air[0]};

    // cerchio che contiene power line
    x_pl = (xx[0]+xx[1]+xx[2])/3;
    y_pl = (yy[0]+yy[1]+yy[2])/3;
    R_pl = Min(25,Abs(y_pl)*0.7); // min(25,10)
    lms_pl = 2.5;
    p_centro = newp; Point(p_centro) = {x_pl, y_pl, 0, lms_pl};
	p_sotto = newp; Point(p_sotto) = {x_pl, y_pl-R_pl, 0, lms_pl};
	p_dx = newp; Point(p_dx) = {x_pl+R_pl, y_pl, 0, lms_pl};
	p_sx = newp; Point(p_sx) = {x_pl-R_pl, y_pl, 0, lms_pl};
    p_sopra = newp; Point(p_sopra) = {x_pl, y_pl+R_pl, 0, lms_pl};

    arc[0] = newl; Circle(arc[0]) = {p_sopra, p_centro, p_sx};
    arc[1] = newl; Circle(arc[1]) = {p_sx, p_centro, p_sotto};
    arc[2] = newl; Circle(arc[2]) = {p_sotto, p_centro, p_dx};
    arc[3] = newl; Circle(arc[3]) = {p_dx, p_centro, p_sopra};

    // SOTTO
    Transfinite Curve {arc[1],arc[2]} = 45 Using Progression 1.0;
    //SOPRA
    Transfinite Curve {arc[3],arc[0]} = 25 Using Progression 1.0;

    ll_int_S_soil[0] = newll;
    Curve Loop(ll_int_S_soil[0]) = {arc[]};

    // SUPERFICIE POWER LINE
    S_pl = news; Plane Surface(S_pl) = {ll_int_S_soil[0]};

	// DIVISIONE TERRENO/ARIA
	linea_terr_dx = newl;
	Line(linea_terr_dx) = {A_ext, p_coppa_p_dx};
	Transfinite Curve {linea_terr_dx} = 200 Using Progression 0.97;
	//
	linea_terr_sx = newl;
	Line(linea_terr_sx) = {p_coppa_p_sx, C_ext};
	Transfinite Curve {linea_terr_sx} = 200 Using Progression 1.03;

	// SUPERF TERRENO/ARIA
	ll_aria = newll;
    // Printf('%g,%g',Cir_B,l_coppa_p_sx);
	Curve Loop(ll_aria) = {Cir_B,linea_terr_sx,l_coppa_p_sx,l_coppa_p_dx,linea_terr_dx,Cir_A};
    S_air = news;
	Plane Surface(S_air) = {ll_aria,ll_int_S_soil[0]}; // ARIA

	ll_terr = newll;
    // Printf('%g,%g',Cir_C,linea_terr_dx);
	Curve Loop(ll_terr) = {Cir_C, Cir_D, linea_terr_dx,l_coppa_p_dx,l_coppa_p_sx,linea_terr_sx};
	S_soil = news;
	Plane Surface(S_soil) = {ll_terr,ll_int_S_air[0]}; // TERRENO

	// PHASE 1
	x = xx[0]; // x centro conduttore
	y = yy[0]; // y centro conduttore
	R_ext = R_ext_fase; // raggio esterno regione strutturata
	R_int = R_int_fase; // raggio interno regione strutturata
	lms = 0.002; // local mesh size pto centrale
	i_physical_surface = 1;
	name_physical_surface = "phase_1";
	flag_isPipe = 0; // se = 1, due physical surf separate per interno/ext conduttore, come pipe
	n_nodi_azimut = 10; // #nodi in direzione azimutale per ogni spicchio di conduttore
	n_nodi_radial = 25; // #nodi in direzione radiale per ogni spicchio di conduttore
	progr_radial = 1.2; // addensamento in direzione radiale per effetto pelle
	// flag_circControllo = 1; // inserisco circonferenza di controllo
	Call GenCondCil;

	// PHASE 2
	x = xx[1]; // x centro conduttore
	y = yy[1]; // y centro conduttore
	i_physical_surface = 2;
	name_physical_surface = "phase_2";
	Call GenCondCil;

	// PHASE 3
	x = xx[2]; // x centro conduttore
	y = yy[2]; // y centro conduttore
	i_physical_surface = 3;
	name_physical_surface = "phase_3";
	Call GenCondCil;

	// OGW
	x = xx[3]; // x centro conduttore
	y = yy[3]; // y centro conduttore
	R_ext = R_ext_OGW;
	R_int = R_int_fase;
	i_physical_surface = 4;
	name_physical_surface = "OGW";
	Call GenCondCil;

	// PIPE
	x = x_pipe_racc(i_spost);
	y = y_pipe;
	R_ext = R_ext_pipe; // raggio esterno regione strutturata
	R_int = R_int_pipe; // raggio interno regione strutturata
	lms = 0.1; // local mesh size pto centrale
	i_physical_surface = 7; // physical: 7(pipe); 8 (centro)
	flag_isPipe = 1; // se = 1, due physical surf separate per interno/ext conduttore, come pipe
	n_nodi_azimut = 30; // #nodi in direzione azimutale per ogni spicchio di conduttore
	n_nodi_radial = 50; // #nodi in direzione radiale per ogni spicchio di conduttore
	progr_radial = 1.2; // addensamento in direzione radiale per effetto pelle
	Call GenCondCil;

	i_mit = 0; // aux
    // L_xx[0] = 3.2; // larghezza lastra
    // L_yy[0] = 0.008; // altezza lastra
	
	// MIT 1R_mit
	x = xx_mit[i_mit]; // x centro conduttore
	y = yy_mit[i_mit]; // y centro conduttore
	R_ext = R_mit; // raggio esterno regione strutturata
	R_int = R_mit/10; // raggio interno regione strutturata
	lms = lms_mit; // local mesh size pto centrale
	i_physical_surface = 9;
	name_physical_surface = "Mit_1";
	flag_isPipe = 0; // se = 1, due physical surf separate per interno/ext conduttore, come pipe
	n_nodi_azimut = 10; // #nodi in direzione azimutale per ogni spicchio di conduttore
	n_nodi_radial = 40; // #nodi in direzione radiale per ogni spicchio di conduttore
	progr_radial = 1.05; // addensamento in direzione radiale per effetto pelle
	Call GenCondCil;
	i_mit ++;
	
	// MIT 2
	x = xx_mit[i_mit]; // x centro conduttore
	y = yy_mit[i_mit]; // y centro conduttore
	i_physical_surface = 10;
	name_physical_surface = "Mit_2";
	Call GenCondCil;
	i_mit ++;

	// MIT 3
	x = xx_mit[i_mit]; // x centro conduttore
	y = yy_mit[i_mit]; // y centro conduttore
	i_physical_surface = 11;
	name_physical_surface = "Mit_3";
	Call GenCondCil;
	i_mit ++;

	// MIT 4
	x = xx_mit[i_mit]; // x centro conduttore
	y = yy_mit[i_mit]; // y centro conduttore
	i_physical_surface = 12;
	name_physical_surface = "Mit_4";
	Call GenCondCil;
	i_mit ++;

	/************************************************************/
	// CIRCONFERENZE DI CONTROLLO
	/************************************************************/
	// FASE 1
	x = xx[0];
	y = yy[0];
	R_ext = R_ext_fase; // raggio esterno regione strutturata
	R_circControllo = 100*R_ext; // raggio circonferenza di controllo
	n_nodi_cirControllo = 25; // # nodi circonferenza di controllo
	Call GenCirControllo;
	R_circControllo = 20*R_ext; // raggio circonferenza di controllo
	Call GenCirControllo;

	// FASE 2
	x = xx[1];
	y = yy[1];
	R_ext = R_ext_fase; // raggio esterno regione strutturata
	R_circControllo = 100*R_ext; // raggio circonferenza di controllo
	n_nodi_cirControllo = 25; // # nodi circonferenza di controllo
	Call GenCirControllo;
	R_circControllo = 20*R_ext; // raggio circonferenza di controllo
	Call GenCirControllo;

	// FASE 3
	x = xx[2];
	y = yy[2];
	R_ext = R_ext_fase; // raggio esterno regione strutturata
	R_circControllo = 100*R_ext; // raggio circonferenza di controllo
	n_nodi_cirControllo = 25; // # nodi circonferenza di controllo
	Call GenCirControllo;
	R_circControllo = 20*R_ext; // raggio circonferenza di controllo
	Call GenCirControllo;

	// OGW
	x = xx[3];
	y = yy[3];
	R_ext = R_ext_OGW; // raggio esterno regione strutturata
	R_circControllo = 150*R_ext; // raggio circonferenza di controllo
	n_nodi_cirControllo = 25; // # nodi circonferenza di controllo
	Call GenCirControllo;
	R_circControllo = 20*R_ext; // raggio circonferenza di controllo
	Call GenCirControllo;

	// PIPE
	x = x_pipe_racc(i_spost);
	y = y_pipe;
	R_ext = R_ext_pipe; // raggio esterno regione strutturata
	R_circControllo = 1.5*R_ext; // raggio circonferenza di controllo
	n_nodi_cirControllo = 100; // # nodi circonferenza di controllo
	Call GenCirControllo;

	i_mit = 0;
	// MIT 1
	x = xx_mit[i_mit];
	y = yy_mit[i_mit];
	L_x = 0.014179631; // base
	L_y = 0.014179631; // altezza
	R_circControllo = 2*Sqrt(L_x^2+L_y^2); // raggio circonferenza di controllo
	n_nodi_cirControllo = 50; // # nodi circonferenza di controllo
	Call GenCirControllo;
	// Call GenEllControllo;
	i_mit ++;
	
	// MIT 2
	x = xx_mit[i_mit];
	y = yy_mit[i_mit];
	Call GenCirControllo;
	i_mit ++;

	// MIT 3
	x = xx_mit[i_mit];
	y = yy_mit[i_mit];
	Call GenCirControllo;
	i_mit ++;

	// MIT 4
	x = xx_mit[i_mit];
	y = yy_mit[i_mit];
	Call GenCirControllo;
	i_mit ++;
	
	/************************************************************/
	
	// ARIA
	i_physical_surface = 5;
	name_physical_surface = "Aria";
	raccolta_superfici = {S_air,S_pl};
	Call Def_physical_surface;
	
	// TERRENO
	i_physical_surface = 6;
	name_physical_surface = "Terreno";
	raccolta_superfici = {S_soil,S_coppa_p};
	Call Def_physical_surface;

	// CONTORNO ARIA 
	b_air() = Boundary{Surface{S_air};};
	Printf("tags AIR: ", b_air()); // stampo tags physycal entities di S_air
	
	i_physical_curve = 71;
	name_physical_curve = "Contorno_aria_1";
	raccolta_curve = {b_air(0)};	
	Call Def_physical_curve;
	
	i_physical_curve = 72;
	name_physical_curve = "Contorno_aria_2";
	raccolta_curve = {b_air(5)};	
	Call Def_physical_curve;

	// CONTORNO SUOLO
	b_soil() = Boundary{Surface{S_soil};};
	Printf("tags SOIL: ", b_soil()); // stampo tags physycal entities di S_soil
	
	i_physical_curve = 73;
	name_physical_curve = "Contorno_suolo_1";
	raccolta_curve = {b_soil(0)};	
	Call Def_physical_curve;
	Printf("tag Contorno_suolo_1: ", b_soil(0));
	
	i_physical_curve = 74;
	name_physical_curve = "Contorno_suolo_2";
	raccolta_curve = {b_soil(1)};	
	Call Def_physical_curve;
	Printf("tag Contorno_suolo_2: ", b_soil(1));
	
	// MESH 2D
	Mesh 2;
	Mesh.MshFileVersion = 2;

	// nome cartella a seconda di posizione pipeline
	If (x_pipe_racc[i_spost] >= 0)
		mydir = Sprintf("MESH_%02g00_pos_%02g",i_spost+1,x_pipe_racc[i_spost]);
	Else
		mydir = Sprintf("MESH_%02g00_neg_%02g",i_spost+1,-x_pipe_racc[i_spost]);
	EndIf

	// creo cartella secondo nome <mydir> e salvo la mesh dentro
	/*
	CreateDir Sprintf(mydir);
	Save Sprintf(StrCat(mydir,"/","mesh.msh"));
	Save Sprintf(StrCat(mydir,"/","mesh.m"));
	*/
	
	Save "mesh_LAPLACE.m";
	
	If (n_spost>1 && i_spost<(n_spost-1) )
		Delete Model;
		Delete Physicals;
	EndIf
	
EndFor