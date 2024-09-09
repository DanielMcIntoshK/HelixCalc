#include <iostream>
#include <string>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <fstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

Matrix cross(Matrix m1, Matrix m2){
	return Matrix{{m1(1,0)*m2(2,0)-m2(1,0)*m1(2,0)},
		{-(m1(0,0)*m2(2,0)-m2(0,0)*m1(2,0))},
		{m1(0,0)*m2(1,0)-m2(0,0)*m1(1,0)}};
}

int main(int argc, char**argv){
	if(argc <5){
		std::cout << "ERROR: Program requires input of distance, angle, dihedral, and number of particels\n";
		return 1;
	}
	
	double d = std::stod(argv[1]);
	double a = std::stod(argv[2]);
	double t = std::stod(argv[3]);
	int n = std::stoi(argv[4]);


	// THIS SECTION JUST GENERATES n POINTS IN A HELIX 
	// IT DOES THIS BY INITIALIZING SOME POINT AT POINT (0,0,0)
	// THEN CREATES A SECOND AT A DISTANCE OF d FROM THE ORIGIN ON THE X AXIS
	// THEN IT CREATES A THIRD AT A DISTANCE OF d FROM THE SECOND WITH AN ANGLE OF a
	// THEN ADDS A FOURTH TAKING INTO ACCOUNT THE TORTIONAL ANGLE
	// IT THEN TRANSLATES ALL THE POINTS SO THAT x2->x1, x3->x2, x4->x3, AND x1 IS PUSHED BACK
	// THEN ANOTHER POINT IS ADDED IN THE SAME SPOT AS THE FOURTH POINT AND THE PATTERN CONTINUES
	double ia =3.141592-a;

	Matrix I=Matrix::Identity(4,4);
	//TRANSLATION MATRIX
	Matrix T{{1.0,0.0,0.0,d},
		{0.0,1.0,0.0,0.0},
		{0.0,0.0,1.0,0.0},
		{0.0,0.0,0.0,1.0}};
	//a ROATION MATRIX
	Matrix A{{std::cos(ia),-std::sin(ia),0.0,0.0},
		{std::sin(ia),std::cos(ia),0.0,0.0},
		{0.0,0.0,1.0,0.0},
		{0.0,0.0,0.0,1.0}};
	//t TORSIONAL ROTATION MATRIX
	Matrix B{{1.0,0.0,0.0,0.0},
		{0.0,std::cos(t),-std::sin(t),0.0},
		{0.0,std::sin(t),std::cos(t),0.0},
		{0.0,0.0,0.0,1.0}};

	//NOTE THE POINTS USE 4 COORDINATES TO ALLOW FOR TRANSLATION MATRICIES TO BE USED
	//THE FORUTH COORDINATE WILL ALWAYS BE 1
	std::cout << d << " " << a << " " << t << " " << n << std::endl;
	Matrix points(4,n);
	points.col(0)=Matrix{{0.0},{0.0},{0.0},{1.0}};
	points.col(1)=T*points.col(0);
	points.col(2)=T*A*T*points.col(0);
	points.col(3)=T*A*B*T*A*T*points.col(0);

	Matrix refpoint=points.col(3);

	for(int i = 4; i<n; i++){
		points=B.inverse()*A.inverse()*T.inverse()*points;
		points.col(i)=refpoint;
	}

	//CALCULATE THE CENTER OF MASS AND TRANSLATE COORDINATES SO IT IS THE ORIGIN
	Matrix sum=Matrix::Zero(4,1);
	for(int i = 0; i < n; i++) sum+=points.col(i);
	sum/=(double)n;
	std::cout << "CENTER OF MASS:\n" <<sum << std::endl;
	
	Matrix cmpoints(3,n);
	for(int i = 0; i < n; i++){
		cmpoints(0,i)=points(0,i)-sum(0,0);
		cmpoints(1,i)=points(1,i)-sum(1,0);
		cmpoints(2,i)=points(2,i)-sum(2,0);
	}

	//CALCULATE THE INTERTIA TENSOR
	Matrix inertiaTensor=Matrix::Identity(3,3);
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			for(int k = 0; k < n; k++){
				if(i==j){
					inertiaTensor(i,i)+=(cmpoints.col(k).transpose()*cmpoints.col(k))(0,0);
				}
				inertiaTensor(i,j)-=cmpoints(i,k)*cmpoints(j,k);
			}
		}
	}
	std::cout <<"INERTIA TENSOR:\n"<< inertiaTensor << std::endl;

	//DIAGONALIZE THE INERTIA TENSOR
	Eigen::SelfAdjointEigenSolver<Matrix> eigen(inertiaTensor);
	std::cout << "THE EIGENVALUES ARE:\n" << eigen.eigenvalues() << std::endl;
	std::cout << "THE EIGENVECTORS ARE:\n" << eigen.eigenvectors() << std::endl;

	Matrix evecs=eigen.eigenvectors();

	//SWAP THE EIGENVECTORS SO LOWEST EIGENVALUE WILL ALLIGN WITH Z AXIS
	int smallest=0;
	for(int i = 1; i< 3;i++){
		if(evecs(i) < evecs(0)) smallest=i;
	}
	if(smallest ==0){
		Matrix ev1=evecs.col(1), ev2=evecs.col(2);
		evecs.col(2)=evecs.col(0);
		evecs.col(1)=ev2;
		evecs.col(0)=ev1;
	}
	if(smallest==1){
		Matrix ev0=evecs.col(0), ev2=evecs.col(2);
		evecs.col(2)=evecs.col(1);
		evecs.col(0)=ev2;
		evecs.col(1)=ev0;
	}
	
	//OUTPUTS THE XYZ FILE OF THE UNTRANSFORMED COORDINATES
	std::ofstream file1("output1.xyz");
	file1 << n << std::endl << "XYZ file made by HelixCalc" << std::endl;
	for(int i = 0; i < n; i++){
		file1 << "H " << cmpoints(0,i) << " " << cmpoints(1,i) << " " << cmpoints(2,i)<<std::endl;
	}
	
	//TRANSLATES COORDINATES TO EIGENBASIS
	cmpoints=evecs.transpose()*cmpoints;

	//OUTPUTS TRANSLATED XYZ FILE
	std::ofstream file2("output2.xyz");
	file2 << n << std::endl << "XYZ file made by HelixCalc" << std::endl;
	for(int i = 0; i < n; i++){
		file2 << "H " << cmpoints(0,i) << " " << cmpoints(1,i) << " " << cmpoints(2,i)<<std::endl;
	}
	
	return 0;
}
