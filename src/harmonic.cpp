#include "harmonic.h"


void Harmonic::star_map() {
	Point center(0, 0, 0);
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		center += vertex->point();
	}
	center /= _nmesh->numVertices();
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
		v_n(vertex) = vertex->point();
	}
}

// e_k(edge) holds cotangent weight.
// edgeVertex1 & edgeVertex2 are edge end points
double Harmonic::compute_energy(int type) {
	double energy = 0;
	for (SolidEdgeIterator eiter(_nmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Solid::tVertex v1 = _nmesh->edgeVertex1(edge);
		Solid::tVertex v2 = _nmesh->edgeVertex2(edge);
		Point uv = v1->point() - v2->point();
		if (type == HARMONIC) {
			energy += e_k(edge) * uv.norm2();
		}
		// Point.norm2() is squared 2-norm
		else energy += uv.norm2();
	}
	return energy;
}

void Harmonic::compute_gradient(int type) {
	// Loop over all vertices sv in mesh
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex sv = *viter;
		// Initialize gradient at sv
		Point gradient = Point(0, 0, 0);
		// Loop over all vertices tv connected with sv
		for (VertexVertexIterator vviter(sv); !vviter.end(); ++vviter) {
			Vertex * tv = *vviter;
			// Gradient computation depends on whether type is HARMONIC or TUETTE:
			// Actually, gradient might be missing a factor 2, given that the cost is ||v1-v2||^2
			if (type == HARMONIC)
				gradient += (sv->point() - tv->point()) * e_k(_nmesh->vertexEdge(sv, tv));
			else gradient += (sv->point() - tv->point());
		}
		// v_dv stores the raw gradient at sv
		v_dv(sv) = gradient;
		// v_abdv stores only the tangential component, since projection onto normal is subtracted
		v_abdv(sv) = v_dv(sv) - v_n(sv) * (v_dv(sv) * v_n(sv));
	}
}

void Harmonic::conjugate_gradient() {
	// ---------------
	// Conclusion: Some doubts about this "conjugate" part, and what it does,
	// but hopefully everything will work out if just modifying compute_gradient.
	// ---------------

	// v_abdv holds tangential component of gradient
	// both v_fabdv AND v_s appear to have been set equal to v_abdv before conjugate gradient loop starts..? Shared value or ref..?
	// Consequently (if shared ref), v_beta might simply be 1, and v_s=0..?
	// v_s redefined to v_s = v_s*v_beta - v_abdv, i.e. incldues both old v_s (weighted), and the negative gradient.
	// In the end, v_s seems to be the direction for vertex update
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		// operator* is dot product, so v_beta(vertex) is set to a squared ratio of vector norms:
		v_beta(vertex) = v_abdv(vertex)*v_abdv(vertex) / (v_fabdv(vertex)*v_fabdv(vertex));
		v_s(vertex) = v_s(vertex)*v_beta(vertex) - v_abdv(vertex);
	}
	double alpha = 1e-6;
	double e0 = compute_energy(HARMONIC);
	// Initialize v_mp to prev value
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v_mp(v) = v->point();
	}
	// Alpha starts small, and is gradually increased, until > 1e-2.
	while (alpha<1e-2) {
		// Perturb the vertices in v_s direction, with step size alpha
		for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
			Solid::tVertex v = *viter;
			v->point() = v->point() + v_s(v)*alpha;
			v->point() /= v->point().norm();
		}
		double nenergy = compute_energy(HARMONIC);
		if (nenergy < e0) {
			// Whenever energy is lower than initial energy, store the value to v_mp. Otherwise, reject.
			for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
				Solid::tVertex v = *viter;
				v_mp(v) = v->point();
			}
		}
		alpha *= 2;
	}
	// v_mp is the final value of the vertex after this conjugate gradient iteration
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v_mp(v);
	}
}

void Harmonic::update_mesh() {
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v->point()- v_abdv(v)*DELTA_T;
		v->point() /= v->point().norm();
	}
}

// Should this be interpreted as an approximate clip operation, in order to fulfill a constraint for sphere being centered at mass center?
void Harmonic::update_mass_center() {
	// Compute per-face areas
	for (SolidFaceIterator fiter(_nmesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;

		Point p[3];
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex * v = *fviter;
			p[i++] = v->point();
		}

		Point n = (p[1] - p[0]) ^ (p[2] - p[0]);
		f_a(face) = n.norm() / 2.0;
	}
	// Compute per-vertex areas (1/3 of surrounding faces)
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		double area = 0;
		for (VertexFaceIterator vfiter(vertex); !vfiter.end(); ++vfiter) {
			Face* vf = *vfiter;
			area += f_a(vf);
		}
		v_a(vertex) = area / 3.0;
	}

	// Increment area-weighted vertex average, and the sum of weights
	Point center(0, 0, 0);
	double mass = 0;
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		center += vertex->point() * v_a(vertex);
		mass += v_a(vertex);
	}
	// Compute the area-weighted vertex average
	center /= mass;

	// Loop over vertices
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		// Translate mass-center to origin
		vertex->point() -= center;
		// Renormalize points
		vertex->point() /= vertex->point().norm();
	}
}

void Harmonic::harmonic_map_conjugate_gd() {
	double te = compute_energy(HARMONIC);
	double oe = 1000;
	compute_gradient(HARMONIC);
	update_mesh();
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		v_fabdv(vertex) = v_abdv(vertex);
		v_s(vertex) = v_abdv(vertex);
	}
	double he = compute_energy(HARMONIC);
	int i = 0;
	while (oe - he > DELTA_E) {
		compute_gradient(HARMONIC);
		conjugate_gradient();
		oe = he;
		he = compute_energy(HARMONIC);
		if (i++ < 10 or i % 100 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::harmonic_map_gd() {
	double he = compute_energy(HARMONIC);
	double oe = 1000;
	int i = 0;
	while (oe - he > DELTA_E) {
		compute_gradient(HARMONIC);
		update_mesh();
		update_mass_center();
		oe = he;
		he = compute_energy(HARMONIC);
		if (i++ < 10 or (i < 500 and i % 100 == 0) or i % 1000 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::tuette_map() {
	double te = compute_energy(TUETTE);
	double oe = 1000;
	int i = 0;
	while (oe-te > DELTA_E) {
		compute_gradient(TUETTE);
		update_mesh();
		oe = te;
		te = compute_energy(TUETTE);
		if (i++ < 10 or (i < 500 and i % 100 == 0) or i % 1000 == 0) {
			std::cout << i << ", energy: " << te << ", diff: "<< oe - te << std::endl;	
		}
	}
}

void Harmonic::set_up_normal() {
	for (SolidVertexIterator viter(_nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->trait() = new CVertexTrait;
	}
	for (SolidFaceIterator fiter(_nmesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;
		face->trait() = new CFaceTrait;
		Point fp[3]; 
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex* fv = *fviter;
			fp[i++] = fv->point();
		}
		Point normal = (fp[1] - fp[0]) ^ (fp[2] - fp[0]);
		f_n(face) = normal / normal.norm();
	}
}

void Harmonic::set_up_kuv() {
	// Initialize "Edge Trait"
	// Dummy initialization of edge weights
	for (SolidEdgeIterator eiter(_nmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		edge->trait() = new CEdgeTrait;
		e_k(edge) = 1.0;
	}

	// Compute and store cotangent weights
	double sum = 0;
	for (SolidEdgeIterator eiter(_nmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Point p1 = _nmesh->edgeVertex1(edge)->point();
		Point p2 = _nmesh->edgeVertex2(edge)->point();
		Point p3 = edge->halfedge(0)->he_next()->target()->point();

		double alpha, beta;
		alpha = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		p3 = edge->halfedge(0)->he_sym()->he_next()->target()->point();
		beta = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		e_k(edge) = alpha + beta;
		sum += e_k(edge);
	}
}

void Harmonic::set_up() {
	set_up_normal();
	set_up_kuv();
}

std::vector<double> Harmonic::compute_vertex_angles() {
	std::vector<double> angles;
	for (SolidFaceIterator fiter(_nmesh); !fiter.end(); ++fiter) {
		Solid::tFace f = *fiter;
		Point p[3];
		int i = 0;
		for (FaceVertexIterator fviter(f); !fviter.end(); ++fviter) {
			Vertex * v = *fviter;
			p[i++] = v->point();
		}
		Point p1 = p[0] - p[1];
		Point p2 = p[1] - p[2];
		Point p3 = p[2] - p[0];

		angles.push_back(acos(p1 * (-p3) / (p1.norm() * (-p3).norm())));
		angles.push_back(acos((-p1) * p2 / ((-p1).norm() * p2.norm())));
		angles.push_back(acos(p3 * (-p2) / (p3.norm() * (-p2).norm())));
	}
	return angles;
}
