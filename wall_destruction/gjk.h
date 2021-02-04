/* This implementation is adapted from an already existing implementation https://github.com/CoppaMan/Pinballs/blob/0db3de619b50ee5fc5f38cccb5bcac7bab1bfc30/8_pinballs/gjk2.h
However, small fixes and additional implementation had to be done
    - Change the vectors from standard vectors to Eigen vectors
    - Adjust input variables
    - handle_thetrahidral did some wrong calculation
    - Only the point vertex point was returned, so vertex to face could return 4 different points
    - Edge to edge was not implemented
    - Face to face was added, but only works if one face lies in the other face.
*/

#ifndef GJK_HEADER
#define GJK_HEADER

#include <Eigen/Core>
#include <igl/per_face_normals.h>
#include <vector>
#include "CollisionDetection.h"


double EPSILON = 1e-4;

class GJK {
public:
    GJK();
    ~GJK();


    static bool run(Eigen::MatrixXd& A, Eigen::MatrixXd& B, std::vector<Contact>& contacts, RigidObject* Abody, RigidObject* Bbody) {
        int GJK_MAXITR = 32;

        std::set<std::tuple<int, int, double, double, double>> contacts_set;

        Eigen::Vector3d new_search_dir = A.row(0).transpose() - B.row(0).transpose();


        Eigen::Vector3d a, b, c, d;
        int num_dim = 0;

        a = support(new_search_dir, A, B); // support in direction of origin
        num_dim++;

        for (int iters = 0; iters < GJK_MAXITR; iters++) {
            // closest point to simplex should make sure that a is free to set it
            // with value of new search direction
            if (closest_point_to_simplex(a, b, c, d, new_search_dir, num_dim)) {
                // collision
                Eigen::Vector3d normal;
                std::vector<Eigen::Vector3d> ps;
                Eigen::Vector3d ea, eb;
                bool type = false; // false ist vertex face, true ist edge edge
                if (EPA(a, b, c, d, A, B, normal, ps, type, ea, eb)) {
                    for (int i = 0; i < ps.size(); i++) {
                        // Here we return add each point found as a contact
                        Contact contact;
                        contact.a = Abody;
                        contact.b = Bbody;
                        contact.p = ps[i];
                        contact.n = normal;
                        contact.ea = ea;
                        contact.eb = eb;
                        contact.type = type ? ContactType::EDGEEDGE : ContactType::VERTEXFACE; 
                        auto ret = contacts_set.insert(std::make_tuple(contact.a->getID(), contact.b->getID(), contact.p[0], contact.p[1], contact.p[2]));
                        if (ret.second) {
                            contacts.push_back(contact);
                        }
                        
                    }
                    return true;
                }
                else {
                    // no collision, which should not be the case for us, as we only look at objects, where collision was detected before
                    return false;
                }
            }
            a = support(new_search_dir, A, B);

            if (a.dot(new_search_dir) < 0) { // -new_search_dir is seperating axis
                return false;
            }
        }

        return false;

    }

    static Eigen::Vector3d handle_line(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d, int& num_dim) {
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d ao = -a;
        Eigen::Vector3d new_search_dir;
        if (ao.dot(ab) > 0) {
            // [a,b]
            num_dim = 3; // origin inbweween ab
            new_search_dir = ab.cross(ao).cross(ab);
            c = b; b = a;

        }
        else {
            // [a]
            num_dim = 2; // origin on side a
            new_search_dir = ao; // in direction of origin
            b = a;
        }
        return new_search_dir;
    }

    static Eigen::Vector3d handle_triangle(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d, int& num_dim) {
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d ac = c - a;
        Eigen::Vector3d ao = -a; // 0 -a
        Eigen::Vector3d abc = ab.cross(ac);

        Eigen::Vector3d new_search_dir;

        if ((abc.cross(ac)).dot(ao) > 0) {
            if (ac.dot(ao) > 0) {
                // [a, c]
                new_search_dir = ac.cross(ao).cross(ac);
                num_dim = 3;
                b = a;
            }
            else {
                //****
                if (ab.dot(ao) > 0) {
                    // [a,b]
                    new_search_dir = ab.cross(ao).cross(ab);
                    c = b; b = a;
                    num_dim = 3;
                }
                else {
                    // [a]
                    new_search_dir = ao;
                    b = a;
                    num_dim = 2;
                }
                // ****
            }
        }
        else {
            if ((ab.cross(abc)).dot(ao) > 0) {
                //****
                if (ab.dot(ao) > 0) {
                    // [a,b]
                    new_search_dir = ab.cross(ao).cross(ab);
                    num_dim = 3;
                    c = b; b = a;
                }
                else {
                    // [a]
                    new_search_dir = ao;
                    num_dim = 2;
                    b = a;
                }
                // ****
            }
            else {
                if (abc.dot(ao) > 0) {
                    // [a,b,c]
                    new_search_dir = abc;
                    num_dim = 4;
                    d = c; c = b; b = a;
                }
                else {
                    // -[a,b,c]
                    new_search_dir = -abc;
                    num_dim = 4;
                    d = c; c = b; b = a;
                }
            }
        }
        return new_search_dir;
    }

    static bool handle_thetrahidral(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d, Eigen::Vector3d& new_search_dir, int& num_dim) {
        /*Eigen::Vector3d ABC = (c-a).cross(b-a); 
        Eigen::Vector3d ACD = (d-a).cross(c-a); 
        Eigen::Vector3d ADB = (b-a).cross(d-a);
        */
        // Here we fixed the bug that the triangles were calculated wrong
        Eigen::Vector3d ABC = (b - a).cross(c - a);
        Eigen::Vector3d ACD = (c - a).cross(d - a);
        Eigen::Vector3d ADB = (d - a).cross(b - a);

        Eigen::Vector3d AO = -a;
        num_dim = 4;

        if (ABC.dot(AO) > 0) {
            d = c; c = b; b = a;
            new_search_dir = ABC;
        }
        else if (ACD.dot(AO) > 0) {
            b = a;
            new_search_dir = ACD;
        }
        else if (ADB.dot(AO) > 0) {
            c = d; d = b; b = a;
            new_search_dir = ADB;
        }
        else {
            // if we reach this point we have enclosed the origin and
            // a collision has occured!!
            return true;
        }

        return false;
    }

    static bool closest_point_to_simplex(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d, Eigen::Vector3d& new_search_dir, int& num_dim) {
        // Johnsons algorithm
        if (num_dim == 1) {
            num_dim = 2;
            b = a;
            new_search_dir = -a;  // 1 point only 1 point can be nearest so we search in direction of origin
            return false;
        }
        if (num_dim == 2) { // line segment
            new_search_dir = handle_line(a, b, c, d, num_dim);
            return false;
        }
        else if (num_dim == 3) { // we're triangle
            new_search_dir = handle_triangle(a, b, c, d, num_dim);
            return false;
        }
        else { // num_dim == 4
            if (handle_thetrahidral(a, b, c, d, new_search_dir, num_dim)) {
                return true; // collision
            }
        }
        return false;
    }


    static Eigen::Vector3d support(const Eigen::Vector3d& v, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
        return support(v, A) - support(-v, B);
    }

    static Eigen::Vector3d support(const Eigen::Vector3d& v, const Eigen::MatrixXd& S) {
        Eigen::VectorXd result = S * v;
        Eigen::MatrixXf::Index   maxIndex;
        result.maxCoeff(&maxIndex);
        return S.row(maxIndex).transpose();
    }

    static Eigen::MatrixXf::Index arg_support(const Eigen::Vector3d& v, const Eigen::MatrixXd& S) {
        Eigen::VectorXd result = S * v;
        Eigen::MatrixXf::Index   maxIndex;
        result.maxCoeff(&maxIndex);
        return maxIndex;
    }

    // Here we added a function, which calculates whether a point is in a face
    static bool inBetween(Eigen::Vector3d& point, std::vector<Eigen::Vector3d>& points) {
        for (int i = 0; i < points.size(); i++) {
            for (int j = i + 1; j < points.size(); j++) {
                Eigen::Vector3d lineDirection = points[j] - points[i];
                double pointProjectOnLine = lineDirection.dot(point) / lineDirection.dot(lineDirection);
                double pointIProjectOnLine = lineDirection.dot(points[i]) / lineDirection.dot(lineDirection);
                double pointJProjectOnLine = lineDirection.dot(points[j]) / lineDirection.dot(lineDirection);

                if (pointProjectOnLine >= pointJProjectOnLine || pointProjectOnLine <= pointIProjectOnLine) {
                    return false;
                }
            }
        }
        return true;
    }

    // Here we added a support function, that returns multiple closest points, if they all have the same distance
    static void support(const Eigen::Vector3d& v, const Eigen::MatrixXd& S, std::vector<Eigen::Vector3d>& points) {
        Eigen::VectorXd result = S * v;
        Eigen::MatrixXf::Index   maxIndex;

        result.maxCoeff(&maxIndex);
        for (int i = 0; i < S.rows(); i++) {
            if (result(i) == result(maxIndex)) {
                points.push_back(S.row(i).transpose());
            }
        }
    }

    // Here we added a function, which calculates the intersection of two lines
    static bool intersect(std::tuple<Eigen::Vector3d, Eigen::Vector3d> listlinesA, std::tuple<Eigen::Vector3d, Eigen::Vector3d> listlinesB, Eigen::Vector3d& ret) {
        Eigen::Vector3d lineAStart = std::get<0>(listlinesA);
        Eigen::Vector3d lineBStart = std::get<0>(listlinesB);
        Eigen::Vector3d lineADirection = std::get<1>(listlinesA) - std::get<0>(listlinesA);
        Eigen::Vector3d lineBDirection = std::get<1>(listlinesB) - std::get<0>(listlinesB);
        if (lineADirection.isApprox(lineBDirection) || lineADirection.isApprox(-lineBDirection) || !std::abs((lineADirection.cross(lineBDirection).dot(lineAStart - lineBStart)) <= 0.0001)) { // margin
            return false;
        }

        Eigen::MatrixXd Atosolve(3, 2);
        Atosolve << lineADirection, lineBDirection;
        Eigen::Vector3d btosolve;
        btosolve = lineBStart - lineAStart;

        Eigen::VectorXd x = Atosolve.colPivHouseholderQr().solve(btosolve);
        ret = lineAStart + lineADirection * x(0);

        return true;
    }


    static bool EPA(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::Vector3d& normal, std::vector<Eigen::Vector3d>& contact_points, bool& type, Eigen::Vector3d& ea, Eigen::Vector3d& eb) {
        const int EPA_MAX_NUM_FACES = 64;
        const int EPA_MAX_NUM_LOOSE_EDGES = 32;
        const int EPA_MAX_NUM_ITERATIONS = 64;

        Eigen::Vector3d faces[EPA_MAX_NUM_FACES][4]; //Array of faces, each with 3 verts and a normal

        //Init with final simplex from GJK

        faces[0][0] = a; faces[0][1] = b; faces[0][2] = c;  faces[0][3] = (b - a).cross(c - a).normalized(); //ABC
        faces[1][0] = a; faces[1][1] = c; faces[1][2] = d;  faces[1][3] = (c - a).cross(d - a).normalized(); //ACD
        faces[2][0] = a; faces[2][1] = d; faces[2][2] = b;  faces[2][3] = (d - a).cross(b - a).normalized(); //ADB
        faces[3][0] = b; faces[3][1] = d; faces[3][2] = c;  faces[3][3] = (d - b).cross(c - b).normalized(); //BDC

        int num_faces = 4;
        int closest_face;

        for (int iterations = 0; iterations < EPA_MAX_NUM_ITERATIONS; iterations++) {
            //Find face that's closest to origin
            float min_dist = faces[0][0].dot(faces[0][3]);
            closest_face = 0;
            for (int i = 1; i < num_faces; i++) {
                float dist = faces[i][0].dot(faces[i][3]);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_face = i;
                }
            }

            //search normal to face that's closest to origin
            Eigen::Vector3d search_dir = faces[closest_face][3];
            Eigen::Vector3d p = support(search_dir, A, B);

            if (p.dot(search_dir) - min_dist < EPSILON) {
                // Convergence (new point is not significantly further from origin)
                // Here usually the one point is returned, we changed that to return multiple points
                normal = -faces[closest_face][3].normalized();
                Eigen::Vector3d mtv = faces[closest_face][3] * p.dot(search_dir); //dot vertex with normal to resolve collision along normal!

                // On which Face is the colliding vertex?
                std::vector<Eigen::Vector3d> pointsA;
                std::vector<Eigen::Vector3d> pointsB;
                support(search_dir, A, pointsA);
                support(-search_dir, B, pointsB);
                Eigen::Vector3d temp;
                // Edge collision if support function on both objects return one edge
                if (pointsA.size() == 2 && pointsB.size() == 2) {
                    if (intersect(std::make_tuple(pointsA[0], pointsA[1]), std::make_tuple(pointsB[0], pointsB[1]), temp)) {
                        contact_points.push_back(temp);
                        ea = pointsA[0] - pointsA[1];
                        eb = pointsB[0] - pointsB[1];
                        type = true; // make it edge edge contact
                        return true;
                    }
                } 
                // If one support function return 4 points, we have a face 
                if (pointsA.size() >= 4) {
                    bool found = false;
                    for (auto point : pointsB) {
                        // check for all points in B if they are in the face of object A
                        if (inBetween(point, pointsA)) {
                            contact_points.push_back(point);
                            found = true;
                        }
                    }
                    if (found) {
                        //normal = -normal;
                        return true;
                    }
                }
                if (pointsB.size() >= 4) {
                    bool found = false;
                    for (auto point : pointsA) {
                        // check for all points in A if they are in th eface of B
                        if (inBetween(point, pointsB)) {
                            contact_points.push_back(point);
                            found = true;
                        }
                    }
                    if (found) {
                        return true;
                    }
                }

                return false;
            }

            Eigen::Vector3d loose_edges[EPA_MAX_NUM_LOOSE_EDGES][2]; //keep track of edges we need to fix after removing faces
            int num_loose_edges = 0;

            //Find all triangles that are facing p
            for (int i = 0; i < num_faces; i++)
            {
                if (faces[i][3].dot(p - faces[i][0]) > 0) //triangle i faces p, remove it
                {
                    //Add removed triangle's edges to loose edge list.
                    //If it's already there, remove it (both triangles it belonged to are gone)
                    for (int j = 0; j < 3; j++) //Three edges per face
                    {
                        Eigen::Vector3d current_edge[2] = { faces[i][j], faces[i][(j + 1) % 3] };
                        bool found_edge = false;
                        for (int k = 0; k < num_loose_edges; k++) //Check if current edge is already in list
                        {
                            if (loose_edges[k][1] == current_edge[0] && loose_edges[k][0] == current_edge[1]) {
                                //Edge is already in the list, remove it
                                //THIS ASSUMES EDGE CAN ONLY BE SHARED BY 2 TRIANGLES (which should be true)
                                //THIS ALSO ASSUMES SHARED EDGE WILL BE REVERSED IN THE TRIANGLES (which
                                //should be true provided every triangle is wound CCW)
                                loose_edges[k][0] = loose_edges[num_loose_edges - 1][0]; //Overwrite current edge
                                loose_edges[k][1] = loose_edges[num_loose_edges - 1][1]; //with last edge in list
                                num_loose_edges--;
                                found_edge = true;
                                k = num_loose_edges; //exit loop because edge can only be shared once
                            }
                        }//endfor loose_edges

                        if (!found_edge) { //add current edge to list
                            // assert(num_loose_edges<EPA_MAX_NUM_LOOSE_EDGES);
                            if (num_loose_edges >= EPA_MAX_NUM_LOOSE_EDGES) break;
                            loose_edges[num_loose_edges][0] = current_edge[0];
                            loose_edges[num_loose_edges][1] = current_edge[1];
                            num_loose_edges++;
                        }
                    }

                    //Remove triangle i from list
                    faces[i][0] = faces[num_faces - 1][0];
                    faces[i][1] = faces[num_faces - 1][1];
                    faces[i][2] = faces[num_faces - 1][2];
                    faces[i][3] = faces[num_faces - 1][3];
                    num_faces--;
                    i--;
                }//endif p can see triangle i
            }//endfor num_faces

            //Reconstruct polytope with p added
            for (int i = 0; i < num_loose_edges; i++)
            {
                // assert(num_faces<EPA_MAX_NUM_FACES);
                if (num_faces >= EPA_MAX_NUM_FACES) break;
                faces[num_faces][0] = loose_edges[i][0];
                faces[num_faces][1] = loose_edges[i][1];
                faces[num_faces][2] = p;
                faces[num_faces][3] = ((loose_edges[i][0] - loose_edges[i][1]).cross(loose_edges[i][0] - p)).normalized();

                //Check for wrong normal to maintain CCW winding
                float bias = 0.000001; //in case dot result is only slightly < 0 (because origin is on face)
                if (faces[num_faces][0].dot(faces[num_faces][3]) + bias < 0) {
                    Eigen::Vector3d temp = faces[num_faces][0];
                    faces[num_faces][0] = faces[num_faces][1];
                    faces[num_faces][1] = temp;
                    faces[num_faces][3] = -faces[num_faces][3];
                }
                num_faces++;
            }
        } //End for iterations
        //printf("EPA did not converge\n");
        //Return most recent closest point
        return false;
    }




};

#endif // GJK_HEADER