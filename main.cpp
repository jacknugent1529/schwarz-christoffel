#include "sc_mappings_disk.hpp"
#include "parameter_selection.hpp"
#include "sc_mappings_class.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

std::pair<std::vector<complex<double>>, std::vector<int>> gen_grid(std::vector<double> B_ks, int theta_divisions, int r_divisions) {
    std::vector<std::complex<double> > zs;
    std::vector<int> line;

    auto delta_r = 1. / ((double) r_divisions);
    auto delta_theta = 2 * M_PI / ((double) theta_divisions);
    std::vector<std::pair<double, int>> thetas;
    for (int j = 0; j < theta_divisions; j++) {
        thetas.push_back({j * delta_theta, 0});
    }
    for (auto theta_over_pi : B_ks) {
        thetas.push_back({theta_over_pi * M_PI, 1});
    }
    std::sort(thetas.begin(), thetas.end());
    const int intervals = 30;

    // evenly spaced lines
    int c = 0;
    for (auto i = 0; i < r_divisions; i++) {
        for (int j = 0; j < thetas.size(); j++) {
            auto [theta, type] = thetas[j];
            delta_theta = (j < thetas.size() - 1) ? (thetas[j+1].first - thetas[j].first) : (2 * M_PI - thetas[j].first);
            if (delta_theta < 1e-5) {continue;} // skip places where uniformly placed thetas overlap with lines to corners
            auto r = i * delta_r;
            for (int k = 0; k <= intervals; k++) {
                double t = ((double) k) / ((double) intervals);
                zs.push_back(std::polar(r + delta_r * t, theta));
                line.push_back(3 * c + type);
                zs.push_back(std::polar(r, theta + delta_theta * t));
                line.push_back(3 * c + 2);
            }
            c++;
        }
    }

    int i = r_divisions;
    for (int j = 0; j < thetas.size(); j++) {
        auto [theta, type] = thetas[j];
        delta_theta = (j < thetas.size() - 1) ? (thetas[j+1].first - thetas[j].first) : (2 * M_PI - thetas[j].first);
        auto r = i * delta_r;
        for (int k = 0; k <= intervals; k++) {
            double t = ((double) k) / ((double) intervals);
            zs.push_back(std::polar(r, theta + delta_theta * t));
            line.push_back(3 * c + 2);
        }
        c++;
    }

    return {zs, line};
}

int main(int argc, char *argv[]) {
    // SOLVING CONSTRAINTS
    // points must be in CW order
    /* std::vector<std::complex<double>> points {{1.,1.}, {1.,-1.}, {-1.,-1.}, {-1.,1.}}; */
    /* std::vector<std::complex<double>> points {{1.,1.}, {2.,-1.}, {-1.,-1.}, {-1.5,-.5}}; */
    /* auto [B_ks, beta_ks, c1, c2] = solve_constraints(points); */

    if (argc < 3) {
        std::cout << "Specify mode of input: either `params` or `points`" << std::endl;
        exit(0);
    }
    std::vector<double> B_ks;
    std::vector<double> beta_ks;
    std::complex<double> c1;
    std::complex<double> c2;
    if (std::string(argv[2]) == "params") {
        int n;
        double x;
        char c;
        std::cin >> n;
        for (int i = 0; i < n; i++) {
            std::cin >> x;
            B_ks.push_back(x);
        }
        for (int i = 0; i < n; i++) {
            std::cin >> x;
            beta_ks.push_back(x);
        }
        double a,b;
        std::cin >> a >> b;
        c1 = {a,b};
        std::cin >> a >> b;
        c2 = {a,b};
    } else if (std::string(argv[2]) == "points") {
        int n;
        double x;
        std::vector<double> as;
        std::vector<complex<double>> points;
        std::cin >> n;
        for (int i = 0; i < n; i++) {
            std::cin >> x;
            as.push_back(x);
        }
        for (int i = 0; i < n; i++) {
            std::cin >> x;
            points.push_back({as[i], x});
        }
        std::tie(B_ks, beta_ks, c1, c2) = solve_constraints(points);
    } else {
        std::cout << "not supported" << std::endl;
        exit(0);
    }
    

    for (auto B : B_ks) {
        std::cout << "B " << B << std::endl;
    }
    for (auto beta : beta_ks) {
        std::cout << "beta " << beta << std::endl;
    }
    std::cout << "c1, c2" << c1 << " " << c2 << std::endl;

    auto ws = setup_ws();

    SchwarzChristoffel<double> sc(B_ks, beta_ks, c1, c2); 
    
    int theta_divisions = 12;
    int r_divisions = 12;
    auto [zs, line] = gen_grid(B_ks, theta_divisions, r_divisions);

    std::string fname = argv[1];
    std::ofstream outfile(fname);

    if (!outfile) {
        std::cout << "Failed to open the file: " << fname << std::endl;
        return 1;
    }

    for (int i = 0; i < zs.size(); i++) {
        auto z = zs[i];
        std::cout << "Current iteration: " << i << "/" << zs.size() << " " << z << "\r" << std::flush;
        auto res = sc(z);
        outfile << z.real() << "," << z.imag() << "," << res.real() << "," << res.imag() << "," << line[i] << std::endl;
    }
}
