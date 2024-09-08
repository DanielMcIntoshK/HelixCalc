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

	double ia =3.141592-a;

	Matrix I=Matrix::Identity(4,4);
	Matrix T{{1.0,0.0,0.0,d},
		{0.0,1.0,0.0,0.0},
		{0.0,0.0,1.0,0.0},
		{0.0,0.0,0.0,1.0}};
	Matrix A{{std::cos(ia),-std::sin(ia),0.0,0.0},
		{std::sin(ia),std::cos(ia),0.0,0.0},
		{0.0,0.0,1.0,0.0},
		{0.0,0.0,0.0,1.0}};
	Matrix B{{1.0,0.0,0.0,0.0},
		{0.0,std::cos(t),-std::sin(t),0.0},
		{0.0,std::sin(t),std::cos(t),0.0},
		{0.0,0.0,0.0,1.0}};


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

	for(int i = 0; i < n; i++){
		points.col(i)=Matrix{{2.5*std::cos(i*3.14159/12)},{2.5*std::sin(i*3.14159/12)},{i*3.0},{1.0}};
	}

	std::cout << "POINTS\n";
	for(int i = 0; i < n; i++){
		std::cout << points.col(i).transpose() << std::endl;
	}
	
	std::cout << "DIFFS\n";
	std::vector<Matrix> diff; diff.resize(n-1);
	for(int i = 0; i < n-1;i++){
		diff[i]=points.col(i+1)-points.col(i);
		std::cout << diff[i].transpose() <<std::endl;
	}
	std::cout << "Distances:\n";
	for(int i = 0; i < n-1; i++){
		double v = (diff[i].transpose()*diff[i])(0,0);
		std::cout << v/d << std::endl;
	}	
	std::cout << "Angles:\n";
	for(int i = 0; i < n-2; i++){
		//std::cout << -diff[i].transpose()*diff[i+1] << std::endl;
		double v = (-diff[i].transpose()*diff[i+1])(0,0);
		std::cout << std::acos(v/(d*d)) << std::endl;
	}
	std::cout << "Dihedrals:\n";
	for(int i = 0; i < n-3; i++){
		Matrix norm1=cross(-diff[0],diff[1]),
		       norm2=cross(-diff[1],diff[2]);
		double length= d*d*std::sin(a);
		double torsion=(norm1.transpose()*norm2)(0,0);
		std::cout << std::acos(torsion/(length*length)) << std::endl;
	}

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
	points=Matrix(0,0);

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

	Eigen::SelfAdjointEigenSolver<Matrix> eigen(inertiaTensor);
	std::cout << "THE EIGENVALUES ARE:\n" << eigen.eigenvalues() << std::endl;
	std::cout << "THE EIGENVECTORS ARE:\n" << eigen.eigenvectors() << std::endl;

	Matrix evecs=eigen.eigenvectors();

	int smallest=0;
	for(int i = 1; i< 3;i++){
		if(eigen.eigenvalues()(i) < eigen.eigenvalues()(0)) smallest=i;
	}
	Matrix smallestvec=evecs.col(smallest);
	std::cout << "SMALLEST VEC\n";
	std::cout << smallestvec.transpose() << std::endl;

	double angle=std::atan2(smallestvec(1,0),smallestvec(0,0));
	Matrix zrot{{std::cos(-angle),-std::sin(-angle),0.0},
		    {std::sin(-angle), std::cos(-angle),0.0},
		    {0.0,	      0.0,            1.0}};
	smallestvec=zrot*smallestvec;
	std::cout << smallestvec.transpose() << std::endl;
	angle=std::atan2(smallestvec(0,0),smallestvec(2,0));
	Matrix yrot{{std::cos(-angle), 0.0,std::sin(-angle)},
		    {0.0,	       1.0,              0.0},
		    {-std::sin(-angle),0.0,std::cos(-angle)}};
	smallestvec=yrot*smallestvec;
	std::cout << smallestvec.transpose() << std::endl;

	std::ofstream file1("output1.xyz");
	file1 << n << std::endl << "XYZ file made by HelixCalc" << std::endl;
	for(int i = 0; i < n; i++){
		file1 << "H " << cmpoints(0,i) << " " << cmpoints(1,i) << " " << cmpoints(2,i)<<std::endl;
	}
	cmpoints=yrot*zrot*cmpoints;
	
	std::ofstream file2("output2.xyz");
	file2 << n << std::endl << "XYZ file made by HelixCalc" << std::endl;
	for(int i = 0; i < n; i++){
		file2 << "H " << cmpoints(0,i) << " " << cmpoints(1,i) << " " << cmpoints(2,i)<<std::endl;
	}
	
	return 0;
}
