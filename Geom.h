#include <array>
#include <cmath>
#include <deque>
#include <random>
#include <set>
#include <unordered_set>
#include <vector>

using namespace std;

[[maybe_unused]] constexpr double EPS = 1E-8;
[[maybe_unused]] constexpr double MAX = 1E9;
[[maybe_unused]] constexpr double RAD_TO_DEG = 180 / M_PI;
[[maybe_unused]] constexpr double DEG_TO_RAD = M_PI / 180;

inline bool eq(double a, double b) {
	return fabs(a - b) < EPS;
}

inline double rnd() {
	return rand() / (RAND_MAX + 1.);
}

class Vec2 {
public:
	double x, y;
	
	Vec2(double x = 0) : x(x), y(x) {}
	
	Vec2(double x, double y) : x(x), y(y) {}
	
	bool operator==(const Vec2& v) const {
		return eq(x, v.x) && eq(y, v.y);
	}
	
	double dot(const Vec2& v) const {
		return x * v.x + y * v.y;
	}
	
	double cross(const Vec2& v) const {
		return x * v.y - y * v.x;
	}
	
	double magn() const {
		return sqrt(x * x + y * y);
	}
	
	double dist(const Vec2& v) const {
		return sqrt((x - v.x) * (x - v.x) + (y - v.y) * (y - v.y));
	}
	
	/* Rotate counterclockwise. */
	Vec2 rotate(double a) const {
		return {x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)};
	}
	
	/* Returns angle range 0 to PI. */
	double angle(const Vec2& v) const {
		return dot(v) / (magn() * v.magn());
	}
};

#define VEC2_OP_V(OP)	Vec2 operator OP(const Vec2& v1, const Vec2& v2) {\
							return {v1.x OP v2.x, v1.y OP v2.y};\
						}

#define VEC2_OP_D(OP)	Vec2 operator OP(const Vec2& v1, double v2) {\
							return {v1.x OP v2, v1.y OP v2};\
						}

VEC2_OP_V(+) VEC2_OP_V(-) VEC2_OP_V(*) VEC2_OP_V(/)

VEC2_OP_D(+) VEC2_OP_D(-) VEC2_OP_D(*) VEC2_OP_D(/)

class Line {
public:
	Vec2 a, b, dir;
	
	Line() = default;
	
	Line(const Vec2& a, const Vec2& b) : a(a), b(b), dir(b - a) {}
	
	bool operator==(const Line& l) const {
		return eq(dir.cross(l.dir), 0) && eq(dir.cross(a - l.a), 0);
	}
	
	double dist(const Vec2& v) const {
		return fabs(dir.cross(v - a)) / dir.magn();
	}
	
	bool parallel(const Line& l) const {
		return eq(dir.cross(l.dir), 0);
	}
	
	Vec2 intersect(const Line& l) const {
		double s1 = (a - l.a).cross(l.dir);
		double s2 = (b - l.a).cross(l.dir);
		return (b * s1 - a * s2) / (s1 - s2);
	}
	
	bool is_intersect(const Line& l) const {
		return fmax(a.x, b.x) - fmin(l.a.x, l.b.x) > -EPS &&
			fmax(l.a.x, l.b.x) - fmin(a.x, b.x) > -EPS &&
			fmax(a.y, b.y) - fmin(l.a.y, l.b.y) > -EPS &&
			fmax(l.a.y, l.b.y) - fmin(a.y, b.y) > -EPS &&
			(l.a - a).cross(dir) * (l.b - a).cross(dir) < EPS &&
			(a - l.a).cross(l.dir) * (b - l.a).cross(l.dir) < EPS;
	}
	
	bool is_strict_intersect(const Line& l) const {
		return fmax(a.x, b.x) - fmin(l.a.x, l.b.x) > EPS &&
			fmax(l.a.x, l.b.x) - fmin(a.x, b.x) > EPS &&
			fmax(a.y, b.y) - fmin(l.a.y, l.b.y) > EPS &&
			fmax(l.a.y, l.b.y) - fmin(a.y, b.y) > EPS &&
			(l.a - a).cross(dir) * (l.b - a).cross(dir) < -EPS &&
			(a - l.a).cross(l.dir) * (b - l.a).cross(l.dir) < -EPS;
	}
	
	bool on_the_line(const Vec2& v) const {
		return eq(dir.cross(v - a), 0);
	}
	
	bool on_the_line_segment(const Vec2& v) const {
		return on_the_line(v) && (v - a).dot(v - b) < EPS;
	}
	
	bool strict_on_the_line_segment(const Vec2& v) const {
		return on_the_line(v) && (v - a).dot(v - b) < -EPS;
	}
	
	bool on_the_left(const Vec2& v) const {
		return dir.cross(v - a) > EPS;
	}
	
	bool on_the_right(const Vec2& v) const {
		return dir.cross(v - a) < -EPS;
	}
	
	bool over_the_line(const Vec2& v) const {
		if (dir.x > EPS) return on_the_left(v);
		if (dir.x < -EPS) return on_the_right(v);
		return false;
	}
	
	bool below_the_line(const Vec2& v) const {
		if (dir.x > EPS) return on_the_right(v);
		if (dir.x < -EPS) return on_the_left(v);
		return false;
	}
};

class Polygon {
public:
	vector<Vec2> vertex;
	
	Polygon() = default;
	
	Polygon(Vec2* v, int s) {
		vertex.insert(vertex.end(), v, v + s);
	}
	
	inline int next(int x) const {
		return (x + 1) % vertex.size();
	}
	
	double perimeter() const {
		double sum = 0;
		size_t size = vertex.size();
		for (int i = 0; i < size; ++i) {
			sum += vertex[i].dist(vertex[next(i)]);
		}
		return sum;
	}
	
	double area() const {
		double sum = 0;
		size_t size = vertex.size();
		for (int i = 0; i < size; ++i) {
			sum += vertex[i].cross(vertex[next(i)]);
		}
		return fabs(sum) / 2;
	}
	
	bool contain(const Vec2& v) const {
		Line line = {v, {MAX, MAX * M_PI}};
		int sum = 0;
		size_t size = vertex.size();
		for (int i = 0; i < size; ++i) {
			Line l = {vertex[i], vertex[next(i)]};
			if (l.on_the_line_segment(vertex[i])) return true;
			sum += line.is_intersect(l);
		}
		return sum % 2;
	}
	
	bool strict_contain(const Vec2& v) const {
		Line line = {v, {MAX, MAX * M_PI}};
		int sum = 0;
		size_t size = vertex.size();
		for (int i = 0; i < size; ++i) {
			Line l = {vertex[i], vertex[next(i)]};
			if (l.on_the_line_segment(vertex[i])) return false;
			sum += line.is_intersect(l);
		}
		return sum % 2;
	}
	
	void compute_convex_hull() {
		int size = (int) vertex.size();
		vector<Vec2> v1 = vertex;
		sort(v1.begin(), v1.end(), [](const Vec2& a, const Vec2& b) {
			return eq(a.x, b.x) ? a.y < b.y : a.x < b.x;
		});
		vector<Vec2*> v2(size + 1);
		int vi = 0;
		for (int i = 0; i < size; ++i) {
			while (vi > 1 && (*v2[vi - 1] - *v2[vi - 2]).cross(v1[i] - *v2[vi - 2]) < EPS) --vi;
			v2[vi++] = &v1[i];
		}
		int k = vi;
		for (int i = size - 2; i >= 0; --i) {
			while (vi > k && (*v2[vi - 1] - *v2[vi - 2]).cross(v1[i] - *v2[vi - 2]) < EPS) --vi;
			v2[vi++] = &v1[i];
		}
		vertex.resize(--vi);
		for (int i = 0; i < vi; ++i) {
			vertex[i] = *v2[i];
		}
	}
};

class DCHVertex {
public:
	Vec2 p;
	static Vec2 c;
	
	DCHVertex(const Vec2& p) : p(p) {}
	
	static double atan(const Vec2& v) {
		return atan2(v.x, v.y);
	}
	
	bool operator<(const DCHVertex& v) const {
	   return atan(p - c) > atan(v.p - c);
	}
};

Vec2 DCHVertex::c;

class DynamicConvexHull {
public:
	using It = multiset<DCHVertex>::iterator;
	
	multiset<DCHVertex> t;
	
	DynamicConvexHull() = default;
	
	DynamicConvexHull(Vec2* v, int s) {
		DCHVertex::c = (v[0] + v[1] + v[2]) / 3;
		for (int i = 0; i < 3; ++i) {
			t.insert(v[i]);
		}
		for (int i = 3; i < s; ++i) {
			insert(v[i]);
		}
	}
	
	It next(const It& it) const {
		It _it = it;
		return ++_it == t.end() ? t.begin() : _it;
	}
	
	It prev(const It& it) const {
		It _it = it;
		return --(_it == t.begin() ? t.end() : _it);
	}
	
	void insert(const Vec2& v) {
		It it = t.insert(v);
		It it_n = next(it);
		It it_p = prev(it);
		if ((it_p->p - it_n->p).cross(v - it_n->p) < EPS) {
			t.erase(it);
			return;
		}
		It past = it_n;
		it_n = next(it_n);
		while ((past->p - it->p).cross(it_n->p - past->p) < EPS) {
			t.erase(past);
			past = it_n;
			it_n = next(it_n);
		}
		past = it_p;
		it_p = prev(it_p);
		while ((past->p - it->p).cross(it_p->p - past->p) > -EPS) {
			t.erase(past);
			past = it_p;
			it_p = prev(it_p);
		}
	}
	
	bool contain(const Vec2& v) const {
		It it = t.lower_bound(v);
		if (it == t.end()) it = t.begin();
		It it_p = prev(it);
		return (it_p->p - it->p).cross(v - it->p) < EPS;
	}
	
	double perimater() const {
		double sum = 0;
		It it = t.begin();
		It past = it;
		++it;
		do {
			sum += it->p.dist(past->p);
			past = it;
			it = next(it);
		} while (past != t.begin());
		return sum;
	}
	
	double area() const {
		double sum = 0;
		It it = t.begin();
		It past = it;
		++it;
		do {
			sum += it->p.cross(past->p);
			past = it;
			it = next(it);
		} while (past != t.begin());
		return fabs(sum) / 2;
	}
};

class PointSet {
public:
	vector<Line> line;
	vector<Vec2> vertex;

	PointSet() = default;
	
	PointSet(Line* l, int s) : line(vector<Line>(l, l + s)) {}
	
	void add(Line* l, int s) {
		line.insert(line.end(), l, l + s);
	}
	
	/* Compute the intersection region of the lines' left sides. */
	void compute_intersect() {
		size_t size = line.size();
		vector<Line> source = line;
		sort(source.begin(), source.end(), [](const Line& a, const Line& b) {
			double atan_a = atan2(a.dir.x, a.dir.y);
			double atan_b = atan2(b.dir.x, b.dir.y);
			if (eq(atan_a, atan_b)) return a.on_the_right(b.a);
			return atan_a < atan_b;
		});
		deque<Line> _line;
		deque<Vec2> _vertex;
		_line.push_back(source[0]);
		for (int i = 1; i < size; ++i) {
			auto& a = source[i];
			auto& b = source[i - 1];
			if (a.parallel(b) && a.dir.dot(b.dir) > 0) continue;
			while (!_vertex.empty() && !a.on_the_left(_vertex.back())) {
				_line.pop_back();
				_vertex.pop_back();
			}
			while (!_vertex.empty() && !a.on_the_left(_vertex.front())) {
				_line.pop_front();
				_vertex.pop_front();
			}
			_vertex.push_back(_line.back().intersect(source[i]));
			_line.push_back(a);
		}
		while (!_vertex.empty() && !_line.front().on_the_left(_vertex.back())) {
			_line.pop_back();
			_vertex.pop_back();
		}
		if (_line.size() > 1) {
			_vertex.push_back(_line.front().intersect(_line.back()));
		}
		vertex.resize(_vertex.size());
		line.resize(_vertex.size());
		copy(_vertex.begin(), _vertex.end(), vertex.begin());
		copy(_line.begin(), _line.end(), line.begin());
	}
	
	void to_polygon(Polygon& p) const {
		p.vertex = vertex;
	}
};

class Circle {
public:
	Vec2 c;
	double r = 0;
	
	Circle() = default;
	
	Circle(const Vec2& c, double r) : c(c), r(r) {}
	
	Circle(const Vec2& v1, const Vec2& v2, const Vec2& v3) {
		double m1 = v1.x * v1.x + v1.y * v1.y;
		double m2 = v2.x * v2.x + v2.y * v2.y;
		double m3 = v3.x * v3.x + v3.y * v3.y;
		double a = (v1.y - v3.y) * (m2 - m1) + (v2.y - v1.y) * (m3 - m1);
		double b = (v3.x - v1.x) * (m2 - m1) + (v1.x - v2.x) * (m3 - m1);
		c = Vec2(a, b) / (2 * v1.cross(v3) + v2.cross(v1) + v3.cross(v2));
		r = c.dist(v1);
	}
	
	Circle(Vec2* v, int s) {
		vector<Vec2> points(v, v + s);
		shuffle(points.begin(), points.end(), default_random_engine(0));
		*this = Circle(points[0], 0);
		for (int i = 1; i < s; ++i) {
			if (c.dist(points[i]) < r + EPS) continue;
			*this = Circle(points[i], 0);
			for (int j = 0; j < i; ++j) {
				if (c.dist(points[j]) < r + EPS) continue;
				c = (points[i] + points[j]) / 2;
				r = c.dist(points[i]);
				for (int k = 0; k < j; ++k) {
					if (c.dist(points[k]) < r + EPS) continue;
					*this = Circle(points[i], points[j], points[k]);
				}
			}
		}
	}
	
	double area() const {
		return M_PI * r * r;
	}
	
	double intersection(const Circle& cir) const {
		double d = c.dist(cir.c);
		if (d >= cir.r + r) return 0;
		if (d <= fabs(cir.r - r)) return r < cir.r ? area() : cir.area();
		double t1 = (r * r + d * d - cir.r * cir.r) / (2 * r * d);
		double t2 = (cir.r * cir.r + d * d - r * r) / (2 * cir.r * d);
		double b = r * sqrt(1 - t1 * t1);
		double a1 = acos(t1) * r * r;
		double a2 = acos(t2) * cir.r * cir.r;
		return a1 + a2 - b * d;
	}
};

class KdNode2 {
public:
	Vec2 v;
	Vec2 min_v;
	Vec2 max_v;
	KdNode2* p = nullptr;
	KdNode2* l = nullptr;
	KdNode2* r = nullptr;
	
	KdNode2(const Vec2& v) : v(v), min_v(v), max_v(v) {}
};

class KdTree2 {
public:
	KdNode2* root = nullptr;
	
	KdTree2(Vec2* v, int s) {
		root = build(v, 0, s, 1, nullptr);
	}
	
	static bool dim_x(const Vec2& a, const Vec2& b) {
		return a.x < b.x;
	}
	
	static bool dim_y(const Vec2& a, const Vec2& b) {
		return a.y < b.y;
	}
	
	static double min_dist(const Vec2& v, const KdNode2& n) {
		bool inside_x = v.x > n.min_v.x && v.x < n.max_v.x;
		bool inside_y = v.y > n.min_v.y && v.y < n.max_v.y;
		if (inside_x && inside_y) return 0;
		if (inside_x) return fmin(v.y - n.max_v.y, n.min_v.y - v.y);
		if (inside_y) return fmin(v.x - n.max_v.x, n.min_v.x - v.x);
		double dist = 0;
		dist = fmin(dist, v.dist(n.min_v));
		dist = fmin(dist, v.dist(n.max_v));
		dist = fmin(dist, v.dist({n.min_v.x, n.max_v.y}));
		dist = fmin(dist, v.dist({n.max_v.x, n.min_v.y}));
		return dist;
	}
	
	static double max_dist(const Vec2& v, const KdNode2& n) {
		double dist = 0;
		dist = fmax(dist, v.dist(n.min_v));
		dist = fmax(dist, v.dist(n.max_v));
		dist = fmax(dist, v.dist({n.min_v.x, n.max_v.y}));
		dist = fmax(dist, v.dist({n.max_v.x, n.min_v.y}));
		return dist;
	}

	KdNode2* build(Vec2* v, int l, int r, bool f, KdNode2* p) {
		if (l >= r) return nullptr;
		int mid = (l + r) / 2;
		nth_element(v + l, v + mid, v + r, f ? dim_x : dim_y);
		KdNode2* n = new KdNode2(v[mid]);
		n->p = p;
		n->l = build(v, l, mid, !f, n);
		n->r = build(v, mid + 1, r, !f, n);
		if (n->l != nullptr) {
			n->min_v.x = fmin(n->min_v.x, n->l->min_v.x);
			n->min_v.y = fmin(n->min_v.y, n->l->min_v.y);
			n->max_v.x = fmax(n->max_v.x, n->l->max_v.x);
			n->max_v.y = fmax(n->max_v.y, n->l->max_v.y);
		}
		if (n->r != nullptr) {
			n->min_v.x = fmin(n->min_v.x, n->r->min_v.x);
			n->min_v.y = fmin(n->min_v.y, n->r->min_v.y);
			n->max_v.x = fmax(n->max_v.x, n->r->max_v.x);
			n->max_v.y = fmax(n->max_v.y, n->r->max_v.y);
		}
		return n;
	}
	
	void farthest(const Vec2& v, KdNode2* n, double& d) const {
		if (n == nullptr || max_dist(v, *n) <= d) return;
		d = fmax(d, v.dist(n->v));
		farthest(v, n->l, d);
		farthest(v, n->r, d);
	}

	void nearest(const Vec2& v, KdNode2* n, double& d) const {
		if (n == nullptr || min_dist(v, *n) >= d) return;
		d = fmin(d, v.dist(n->v));
		nearest(v, n->l, d);
		nearest(v, n->r, d);
	}
	
	void farthest_pair(KdNode2* n, double& d) const {
		if (n == nullptr) return;
		farthest(n->v, root, d);
		farthest_pair(n->l, d);
		farthest_pair(n->r, d);
	}
	
	void nearest_pair(KdNode2* n, double& d) const {
		if (n == nullptr) return;
		Vec2 t = n->v;
		n->v.x += d * 2;
		nearest(t, root, d);
		n->v = t;
		nearest_pair(n->l, d);
		nearest_pair(n->r, d);
	}
};

class Vec3 {
public:
	double x, y, z;
	
	Vec3(double x = 0) : x(x), y(x), z(x) {}
	
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	
	bool operator==(const Vec3& v) const {
		return eq(x, v.x) && eq(y, v.y) && eq(z, v.z);
	}
	
	double dot(const Vec3& v) const {
		return x * v.x + y * v.y + z * v.z;
	}
	
	Vec3 cross(const Vec3& v) const {
		return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
	}
	
	double magn() const {
		return sqrt(x * x + y * y + z * z);
	}
	
	double dist(const Vec3& v) const {
		return sqrt((x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z));
	}
	
	Vec3 rotate(const Vec3& v, double a) const {
		return {(cos(a) + (1 - cos(a)) * v.x * v.x) * x +
			((1 - cos(a)) * v.x * v.y - sin(a) * v.z) * y +
			((1 - cos(a)) * v.x * v.z + sin(a) * v.y) * z,
			((1 - cos(a)) * v.x * v.y + sin(a) * v.z) * x +
			(cos(a) + (1 - cos(a)) * v.y * v.y) * y +
			((1 - cos(a)) * v.y * v.z - sin(a) * v.x) * z,
			((1 - cos(a)) * v.x * v.z - sin(a) * v.y) * x +
			((1 - cos(a)) * v.y * v.z + sin(a) * v.x) * y +
			(cos(a) + (1 - cos(a)) * v.z * v.z) * z};
	}
	
	double angle(const Vec3& v) const {
		return dot(v) / (magn() * v.magn());
	}
};

#define VEC3_OP_V(OP)	Vec3 operator OP(const Vec3& v1, const Vec3& v2) {\
							return {v1.x OP v2.x, v1.y OP v2.y, v1.z OP v2.z};\
						}

#define VEC3_OP_D(OP)	Vec3 operator OP(const Vec3& v1, double v2) {\
							return {v1.x OP v2, v1.y OP v2, v1.z OP v2};\
						}

VEC3_OP_V(+) VEC3_OP_V(-) VEC3_OP_V(*) VEC3_OP_V(/)

VEC3_OP_D(+) VEC3_OP_D(-) VEC3_OP_D(*) VEC3_OP_D(/)

class ConvexHull {
public:
	std::vector<Vec3> vertices;
	std::vector<Vec3> normals;
	std::vector<std::array<int, 3> > faces;
	
	ConvexHull(Vec3* v, int s) : vertices(vector<Vec3>(v, v + s)) {}
	
	void add_vertex(Vec3* v, int s) {
		vertices.insert(vertices.end(), v, v + s);
	}
	
	void compute_convex_hull() {
		insert_face(0, 1, 2);
		insert_face(2, 1, 0);
		size_t size = vertices.size();
		for (int i = 3; i < size; ++i) {
			std::unordered_set<long long> new_faces;
			auto face_iter = faces.begin();
			auto normal_iter = normals.begin();
			while (face_iter != faces.end()) {
				Vec3 dir = vertices[i] - vertices[(*face_iter)[0]];
				if (normal_iter->dot(dir) <= 0) {
					++face_iter;
					++normal_iter;
					continue;
				}
				for (int k = 0; k < 3; ++k) {
					long long u = (*face_iter)[k];
					long long v = (*face_iter)[(k + 1) % 3];
					long long uv = v << 32 | u;
					if (new_faces.count(uv) != 0) {
					   new_faces.erase(uv);
					} else {
					   new_faces.insert(u << 32 | v);
					}
				}
				face_iter = faces.erase(face_iter);
				normal_iter = normals.erase(normal_iter);
			}
			for (auto& f : new_faces) {
				insert_face(i, f >> 32, f & 0xFFFFFFFFll);
			}
		}
	}
	
	void insert_face(int a, int b, int c) {
		faces.push_back({a, b, c});
		normals.push_back((vertices[b] - vertices[a]).cross(vertices[c] - vertices[a]));
	}
};

class Ray {
public:
	Vec3 o, dir;
	
	Ray() = default;
	
	Ray(const Vec3& o, const Vec3& d) : o(o), dir(d) {}
	
	/* Returns the distance between o and triangle. */
	double intersect_triangle(const Vec3& a, const Vec3& b, const Vec3& c, Vec3* i = nullptr) {
		Vec3 ab = b - a;
		Vec3 ac = c - a;
		Vec3 ao = o - a;
		Vec3 p = dir.cross(ac);
		Vec3 q = ao.cross(ab);
		double inverse = 1 / ab.dot(p);
		if (i != nullptr) {
			i->y = ao.dot(p) * inverse;
			i->z = dir.dot(q) * inverse;
			i->x = 1 - i->y - i->z;
		}
		return ac.dot(q) * inverse;
	}
};

double simpson(double l, double r, double f(double)) {
	return (f(l) + f(r) + f((l + r) / 2) * 4) * (r - l) / 6;
}

double adaptive_simpson(double l, double r, double s, double p, double f(double)) {
	double mid = (l + r) / 2;
	double sl = simpson(l, mid, f);
	double sr = simpson(mid, r, f);
	if (fabs(sl + sr - s) < p) return sl + sr;
	return adaptive_simpson(l, mid, sl, p, f) + adaptive_simpson(mid, r, sr, p, f);
}
