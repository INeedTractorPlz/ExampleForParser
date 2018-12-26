#include<MyLib/basic_types.hpp>
//#include<MyLib/functions.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include<vector>

#include<fstream>
#include<iostream>
#include<ios>


using namespace boost::numeric::ublas;

typedef double data_type;
typedef unsigned uint;
typedef vector<data_type> vector_data;
typedef matrix<data_type> matrix_data;

vector_data cross_product(const vector_data &u, const vector_data &v){
    vector_data result(u.size());
    result(0) = u(1)*v(2) - u(2)*v(1);
    result(1) = u(2)*v(0) - u(0)*v(2);
    result(2) = u(0)*v(1) - u(1)*v(0);
    return result;
}

struct Newton_t{
    std::vector<data_type>& m;
    uint dim;
    data_type G, pi=atan(data_type(1.))*4;
    uint number_bodies = m.size();
    Newton_t(std::vector<data_type> &m, uint dim, data_type G) : m(m), dim(dim), G(G) { }
    void operator()(const matrix_data& R, matrix_data& A,  data_type t){
        A=zero_matrix<data_type>(number_bodies,dim);
        subrange(A,0,number_bodies,0,dim/2)=subrange(R,0,number_bodies,dim/2,dim);
        matrix_data norm_r(number_bodies,number_bodies);

        norm(R,norm_r);        
        for(uint j=0;j<number_bodies;++j){
            for(uint k=0;k<j;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
            for(uint k=j+1;k<number_bodies;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
        }
        subrange(A,0,number_bodies,dim/2,dim)*=G;
    }    
    
    void norm(const matrix_data& R, matrix_data& norm_r){
        for(uint i=0;i<R.size1();++i){
            for(uint j=0;j<i;++j)
                norm_r(i,j)=norm_r(j,i);
            for(uint j=i+1;j<R.size1();++j){
                norm_r(i,j)=norm_2(row(subrange(R,0,number_bodies,0,dim/2),i)
                -row(subrange(R,0,number_bodies,0,dim/2),j));
            }
        }
    }
    void get_orbital_elements(const vector_data &R, uint i, uint j,
        vector_data& orbital_elements,const matrix_data &norm_r){
        data_type kappa_quad,r,v;
        data_type a,e,inclination;

        
        kappa_quad=G*(m[i]+m[j]);
        r=norm_r(i,j);
        v=norm_2(subrange(R,dim/2,dim));
        
        data_type h = v*v/2 - kappa_quad/r;
        
        vector_data c(3);
        if(dim == 6)
            c = cross_product(subrange(R,0,dim/2),subrange(R,dim/2,dim));
        else{
            c(0) = c(1) = 0;
            c(2) = R(0)*R(3) - R(1)*R(2); 
        }
        
        if(2.*h*inner_prod(c,c)/(kappa_quad*kappa_quad) + 1 < 0){
            assert(2.*h*inner_prod(c,c)/(kappa_quad*kappa_quad) + 1 > -1.0e-10 && "Negative sqrt");
            e = 0;
        }else
            e = sqrt(1.+2.*h*inner_prod(c,c)/(kappa_quad*kappa_quad));
        a = -kappa_quad/2/h;
    
        if(norm_2(c) != 0) 
            inclination = acos(c(2)/norm_2(c));
        else
            inclination = 0;
        
        orbital_elements(0) = a;
        orbital_elements(1) = e;
        orbital_elements(2) = inclination*180/pi;
        orbital_elements(3) = (1-e)*a;
    }
};

struct write_elements
{
    std::ofstream &elements_file;
    std::vector<data_type> &M;
    vector_data elements;
    uint freq;
    uint counter = 0;
    write_elements(std::ofstream &elements_file, std::vector<data_type> &M, uint freq, 
    vector_data elements)
    : elements_file(elements_file), M(M), freq(freq), elements(elements) {}
    void operator()(const state_matrix &x, data_type t)
    {
        if (counter % freq == 0)
            {
                elements_file << std::setw(20) << t << std::setw(20) << M[1];
                for (size_t k = 0; k < elements.size(); ++k)
                    elements_file << std::setw(20) << elements[k];
                    elements_file << std::endl;
                }
                counter++;
    }
};

struct Integrator_t{
    Newton_t Force;
    write_elements &Observer;
    data_type dM = 1.0e-9, dt;
    uint i, j;
    Integrator_t(write_elements &Observer, Newton_t Force, uint i, uint j,
    data_type dt) : Observer(Observer), Force(Force), i(i), j(j), dt(dt) { }
    
    void operator()(const state_matrix& R, state_matrix& A,  data_type t){
        Force(R, A, t);
    }
    void operator()(const state_matrix &R, data_type t){
        matrix_data norm_r(R.size1(), R.size1());
        Force.norm(R, norm_r);
        //simple_cout(R);
        Force.get_orbital_elements(row(R, i) - row(R, j), i, j, Observer.elements, norm_r);
        //simple_cout("Observer.elements = ", Observer.elements);
        //simple_cout("M = ", Observer.M);
        Observer(R, t);
        //simple_cout("after_Obsrever_Yes");
        Force.m[1] -= Force.m[1] * dM * dt;
    }

};

int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    std::vector<data_type> m(2);
    uint dim = 4, number_steps;
    data_type dt, pi=atan(data_type(1.))*4, G=4*pi*pi, T;
    
    std::ofstream elements_file;
    std::ifstream initial_parameters, initial_data;
    uint freq;

    elements_file.precision(12);
    elements_file.open("Result_Eskin_eu.dat", std::ios_base::trunc);
    matrix_data initial(2,dim);

    initial_parameters.open("Parameters_Eskin.dat", std::ios_base::in);
    initial_parameters >> m[0] >> m[1] >> dim >> freq >> T >> number_steps;
    initial_parameters.close();

    initial_data.open("Data_Eskin.dat", std::ios_base::in);
    initial_data >> initial(0,0) >> initial(0,1) >> initial(0,2) >> initial(0,3);
    initial_data >> initial(1,0) >> initial(1,1) >> initial(1,2) >> initial(1,3);
    initial_data.close();

    dt = T/number_steps;

    std::vector<struct option> longOpts({
        {"number_steps",required_argument, NULL,0}
    });

    std::string short_opts = "0:1:d:h:f:T";
    Parser_t<6,1> Parser(short_opts, &longOpts[0]);

    Parser.Parser(argc, argv, m[0], m[1], dim, dt, freq, T, number_steps);

    
    simple_cout("m = ", m);
    simple_cout(sqrt(G*(m[0]+m[1])/40));

    write_elements Observer(elements_file, m, freq, vector_data(4));
    Newton_t Force(m , dim, G);
    //simple_cout("Force.m = ", Force.m);
    Integrator_t Integrator(Observer, Force, 1, 0, dt);
    //simple_cout("Integrator.Froce.m = ", Integrator.Force.m);

    RungeKutta5_Fehlberg<data_type, matrix_data> rk5;
    RungeKutta4<data_type, matrix_data> rk4;
    number_steps = integrate(rk5, std::ref(Integrator), initial, 0.0, T/number_steps, number_steps,
    std::ref(Integrator));

    simple_cout("number_steps = ", number_steps);
    elements_file.close();
    return 0;
}