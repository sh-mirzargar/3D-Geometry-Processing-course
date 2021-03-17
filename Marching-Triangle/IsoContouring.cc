#include "IsoContouring.hh"
#include <math.h>

void
IsoContouring::show_isovalue_and_level_set() {
    show_isovalue();
    create_level_set0_segments();
}

//====================================================================================================================//

void
IsoContouring::show_isovalue() {
    std::vector<double> iso_values;
    for(auto vh : trimesh_.vertices()) {
        auto point = trimesh_.point(vh);
        auto iv = iso_value(point);

        iso_values.push_back(iv);
    }

    auto max_iv = *std::max_element(iso_values.begin(), iso_values.end());
    auto min_iv = *std::min_element(iso_values.begin(), iso_values.end());

    const double range = max_iv - min_iv;
    color_.set_range(0, 1.0, false);

for (auto vh : trimesh_.vertices()) {
    auto t = (iso_values[vh.idx()] - min_iv) / range;
    trimesh_.set_color(vh, color_.color_float4(t));
}
}


double
IsoContouring::iso_value(const Point& _pt) const {
    auto x = _pt[0];
    auto y = _pt[1];
    double iso = 0.;

    //Four functions are given in the handout
    //----------Add your code here----------
    switch (function_id_) {
    case 0: {
        //Function 1
        iso = sqrt((pow(x, 2) + pow(y, 2))) - 1;
        break;
    }
    case 1: {
        //Function 2 f(x,y) = y**2 - sin (x**2) = 0
        iso = pow(y, 2) - sin(pow(x, 2));
        break;
    }
    case 2: {
        //Function 3 f(x,y) = sin(2x+2y)-cos(4xy)+1 = 0
        iso = sin(2 * x + 2 * y) - cos(4 * x * y) + 1;
        break;
    }
    default: {
        //Function 4 f(x,y) = (3*x**2-y**2)**2*y**2-(x**2+y**2)=0
        iso = (pow(3 * pow(x, 2) - pow(y, 2), 2) * pow(y, 2)) - pow((pow(x, 2) + pow(y, 2)), 4);
        break;
    }
    }
    //--------------------------------------

    return iso;
}


void
IsoContouring::create_level_set0_segments() {
    compute_segment_points();

    //create segments
    std::vector<OpenVolumeMesh::VertexHandle> vhs;
    for (auto pt : segment_points_)
        vhs.push_back(line_mesh_.add_vertex(pt));

    if (vhs.size() % 2 != 0)
        std::cout << "\nError: one segment has only one point!";

    for (size_t i = 0; i < vhs.size() / 2; ++i)
        line_mesh_.add_edge(vhs[2 * i], vhs[2 * i + 1]);
}
int
IsoContouring::find_triangle_type(std::vector<double> iso_vec, int type) {
    //int type = 1;
    //std::cerr << iso_vec[0];
    //std::cerr << iso_vec[1];
    //std::cerr << iso_vec[2];
    if (iso_vec[0] > 0) {
        if (iso_vec[1] > 0) {
            if (iso_vec[2] > 0) {
                type = 1; //all positive
            }
            else if (iso_vec[2] < 0) {
                type = 2; //iso[0,1] positive and iso[2]negative
            }
        }
        else if (iso_vec[1] < 0) {
            if (iso_vec[2] > 0) {
                type = 3; //iso[0,2] posiitive and iso[1] negative
            }
            else if (iso_vec[2] < 0) {
                type = 4; //iso[0] positive and iso[1,2] negative
            }
        }
    }
    else if (iso_vec[0] < 0) {
        if (iso_vec[1] > 0) {
            if (iso_vec[2] > 0) {
                type = 5; //iso[0] negative and iso[1,2] positive
            }
            else if (iso_vec[2] < 0) {
                type = 6; //iso[0,2] negative and iso[1] positive
            }
        }
        else if (iso_vec[1] < 0) {
            if (iso_vec[2] > 0) {
                type = 7; //iso [0,1] negative and iso[3] positive
            }
            else if (iso_vec[2] < 0) {
                type = 8; // all negative
            }
        }
    }
    return type;
}

void
IsoContouring::compute_segment_points() {
    segment_points_.clear();
    Point seg_0;
    Point seg_1;
    std::vector<Point> v_points(trimesh_.n_vertices());
    std::vector<std::vector<int> > triangle_ids;

    //store vertex position in v_points
    for (auto vh : trimesh_.vertices()) {
        Point v_pos = trimesh_.point(vh);
        v_points[vh.idx()] = v_pos;
    }

    //store vertex indices in triangle_ids
    for (auto fh : trimesh_.faces()) {
        std::vector<int> vv(3);
        int k = 0;
        for (auto fv_it = trimesh_.fv_ccwbegin(fh); fv_it != trimesh_.fv_ccwend(fh); ++fv_it) {
            vv[k] = (*fv_it).idx();
            k++;
        }
        triangle_ids.push_back(vv);
    }

    
    
    //----------Add your code here----------
    //For each vertex you should compute the scalar iso-value of an implicit function
    //iterate over the triangle indices
    for (std::vector<std::vector<int>>::iterator it = triangle_ids.begin(); it != triangle_ids.end(); ++it) {
        
        
        std::vector<int> triangle_vertexes = *it;
        //for each triangle fetch its vertices
        std::vector<double> iso_vec;
        std::vector<Point> virtices;
        // Then for each triangle, store the iso-value of its virtices and vertex points in two vectors
        for (std::vector<int>::iterator it_2 = triangle_vertexes.begin(); it_2 != triangle_vertexes.end(); ++it_2) {
            iso_vec.push_back(iso_value(v_points[*it_2]));
            virtices.push_back(v_points[*it_2]);
            
        }
        //find triangle_type_fucntion returns 1 if all the vertices have positive iso values,
        //returns 8 if all virtices have negative value and etc.
       
        int type = 0;
        type = find_triangle_type(iso_vec,type);
        //std::cerr << type;
        //If the signs are different, use linear interpolation based on the iso-values 
        //to compute the edge that passes through that triangle   
        
        if (type != 1 && type != 8) {
            
            //coordinates of the point whose sign is unique
            double x_0 = 0;
            double y_0 = 0;
            //coordinates of the points whose sign are equall
            double x_1 = 0;
            double y_1 = 0;

            double x_2 = 0;
            double y_2 = 0;

            double f_0 = 0;
            double f_1 = 0;
            double f_2 = 0;

            //deponding on the type, initialize the x_{0,1,2} and y_{0,1,2}
            switch (type) {
            case 2:
                //iso[0,1] positive and iso[2]negative
                x_0 = v_points[2][0];
                y_0 = v_points[2][1];

                x_1 = v_points[0][0];
                y_1 = v_points[0][1];

                x_2 = v_points[1][0];
                y_2 = v_points[1][1];

                f_0 = iso_vec[2];
                f_1 = iso_vec[0];
                f_2 = iso_vec[1];
                break;
            case 3:
                //iso[0,2] posiitive and iso[1] negative
                x_0 = v_points[1][0];
                y_0 = v_points[1][1];

                x_1 = v_points[0][0];
                y_1 = v_points[0][1];

                x_2 = v_points[2][0];
                y_2 = v_points[2][1];

                f_0 = iso_vec[1];
                f_1 = iso_vec[0];
                f_2 = iso_vec[2];
                break;
            case 4:
                //iso[0] positive and iso[1,2] negative
                x_0 = v_points[0][0];
                y_0 = v_points[0][1];

                x_1 = v_points[1][0];
                y_1 = v_points[1][1];

                x_2 = v_points[2][0];
                y_2 = v_points[2][1];

                f_0 = iso_vec[0];
                f_1 = iso_vec[1];
                f_2 = iso_vec[2];
                break;
            case 5:
                //iso[0] negative and iso[1,2] positive
                x_0 = v_points[0][0];
                y_0 = v_points[0][1];

                x_1 = v_points[1][0];
                y_1 = v_points[1][1];

                x_2 = v_points[2][0];
                y_2 = v_points[2][1];

                f_0 = iso_vec[0];
                f_1 = iso_vec[1];
                f_2 = iso_vec[2];
                break;
            case 6:
                //iso[0,2] negative and iso[1] positive
                x_0 = v_points[1][0];
                y_0 = v_points[1][1];

                x_1 = v_points[0][0];
                y_1 = v_points[0][1];

                x_2 = v_points[2][0];
                y_2 = v_points[2][1];

                f_0 = iso_vec[1];
                f_1 = iso_vec[0];
                f_2 = iso_vec[2];
                break;
            case 7:
                //iso [0,1] negative and iso[2] positive
                x_0 = v_points[2][0];
                y_0 = v_points[2][1];

                x_1 = v_points[0][0];
                y_1 = v_points[0][1];

                x_2 = v_points[1][0];
                y_2 = v_points[1][1];

                f_0 = iso_vec[2];
                f_1 = iso_vec[0];
                f_2 = iso_vec[1];
                break;
            default:
                std::cerr << type <<std::endl;
                std::cerr << "incorrect type number" << std::endl;
            }
            //find the edge that passes through the triangle using linear intrapolation
            //Linear intrapolation implementation is from this site (slide 20):
           //https://www.cse.wustl.edu/~taoju/cse554/lectures/lect04_Contouring_I.pdf
            
            double x_seg_0 = x_0 + (f_0 / (f_0 - f_1) * (x_1 - x_0));
            double y_seg_0 = y_0 + (f_0 / (f_0 - f_1) * (y_1 - y_0));

            double x_seg_1 = x_0 + (f_0 / (f_0 - f_2) * (x_2 - x_0));
            double y_seg_1 = y_0 + (f_0 / (f_0 - f_2) * (y_2 - y_0));

            
            seg_0[0] = x_seg_0;
            seg_0[1] = y_seg_0;
            seg_0[2] = 0;

            seg_1[0] = x_seg_1;
            seg_1[1] = y_seg_1;
            seg_1[2] = 0;

            //segment_points_ is defined in IsoContouring.hh as std::vector<Point> segment_points_;
            //add points in segment_points forming an edge one after the other;
            //for example segment_points[0] and segment_points[1] are two points forming the first edge
            //and segment_points[2] and segment_points[3] are two points forming the second edge
            segment_points_.push_back(seg_0);
            segment_points_.push_back(seg_1);
            std::cerr << "seg" << std::endl;
            std::cerr << seg_0 << std::endl;
            std::cerr << seg_1 << std::endl;
        }
    
           
             
    }
   
    //--------------------------------------

}
