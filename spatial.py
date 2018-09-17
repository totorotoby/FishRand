from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString


class HotSpot:

    def __init__(self, name, deftype, assos_fish, attraction, list):

        self.name = name
        self.deftype = deftype
        self.fish = assos_fish
        self.attraction = attraction
        self.defintion = list
        self.parse_list()
        self.polygon = None
        self.weights = []

    def parse_list(self):

        if self.deftype == 'Polygon':

            for i in range(len(self.defintion)):
                self.defintion[i] = [float(j) for j in self.defintion[i].replace(' ', '').split(',')]

    def addpoly(self, polygon):

        self.polygon = polygon

    def calcweights(self, reg_polys):

        for poly in reg_polys:
            area_of_intersect = poly[1].intersection(self.polygon).area
            self.weights.append([poly[0], area_of_intersect])

        weight_sum = sum(row[1] for row in self.weights)

        for i in range(len(self.weights)):
            unnormed_prob = self.weights[i][1]
            self.weights[i][1] = unnormed_prob/weight_sum

    def get_chem_region(self):

        names = [row[0] for row in self.weights]
        prob = [row[1] for row in self.weights]

        return np.random.choice(names, p=prob)



# From Sergii Nechuiviter https://gist.github.com/Sklavit/e05f0b61cb12ac781c93442fbea4fb55
def voronoi_finite_polygons_2d(vor, radius=None):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

# prints the intersection boundary and regional polygons of input
def print_intersections(polygon, boundary, intersection = None):

    x_b, y_b = boundary.exterior.xy
    if type(polygon) == Polygon:
        x_p, y_p = polygon.exterior.xy
    else:
        coords = polygon.boundary
        x_p = [row.x for row in coords]
        y_p = [row.y for row in coords]

    plt.plot(x_b, y_b, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    plt.plot(x_p, y_p, color='r', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    if intersection != None:
        x_i, y_i = intersection.exterior.xy
        plt.plot(x_i, y_i, color='g', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)

    plt.show()

# Plots a single region
def plot_shape(shapes):

    min_x =[]
    max_x = []
    min_y = []
    max_y = []
    for shape in shapes:

        x_s, y_s = shape.exterior.xy

        x_max = max(x_s)
        bit = x_max * .1
        x_max = x_max + bit
        max_x.append(x_max)
        x_min = min(x_s) - bit
        min_x.append(x_min)
        y_max = max(y_s) + bit
        max_y.append(y_max)
        y_min = min(y_s) - bit
        min_y.append(y_min)

        plt.plot(x_s, y_s, color='black', linewidth=3, solid_capstyle='round', zorder=2)
        plt.xlim((min(min_x), max(max_x)))
        plt.ylim(((min(min_y), max(max_y))))

    plt.show()


def plot_vor(vor, boundary, hotspots):

    fig = voronoi_plot_2d(vor)
    x_b, y_b = boundary.exterior.xy
    x_max = max(x_b)
    bit = x_max * .1
    x_max = x_max + bit
    x_min = min(x_b) - bit
    y_max = max(y_b) + bit
    y_min = min(y_b) - bit
    plt.plot(x_b, y_b, color='black', linewidth=3, solid_capstyle='round', zorder=2)
    plt.xlim((x_min, x_max))
    plt.ylim(((y_min, y_max)))

    for hotspot in hotspots:
        x_h, y_h = hotspot.polygon.exterior.xy
        plt.plot(x_h, y_h, color='g', linewidth=3, solid_capstyle='round', zorder=2)
    plt.show()


# Returns polygons representing the boundary regions and hotspots
def setup(site_data):

    # get boundary polygon
    boundary = Polygon(site_data[0])

    # get regional polygons...takes some messing around with Voronoi function...
    points = np.array([row[1] for row in site_data[1]])
    vor = Voronoi(points)
    regions, vertices = voronoi_finite_polygons_2d(vor, radius=boundary.length * 2)
    reg_polygons = []
    for i in range(len(regions)):
        if len(regions[i]) != 0:
            poly_u = Polygon(vertices[regions[i]])
            poly = poly_u.intersection(boundary)

            # use to see intersections
            #print_intersections(poly_u,boundary, poly)

            reg_polygons.append([site_data[1][i][0],poly])


    attract_poly = site_data[2]
    for hotspot in attract_poly:
        if hotspot.deftype == 'Polygon':
            polygon = Polygon(hotspot.defintion)
            hotspot.addpoly(polygon)
            hotspot.calcweights(reg_polygons)


    # Use to see regions
    # plot_polygons(vor,reg_polygons,vertices, points)
    # plot_vor(vor, boundary, attract_poly)

    return boundary, reg_polygons, attract_poly

# returns probability of being in each hotspot, and probability of being outside of hotspots
def hotspot_prob(boundary, attraction_polys):

    probs = []

    outside = boundary
    for poly in attraction_polys:
        outside = outside.difference(poly.polygon)
        probs.append(poly.polygon.area * poly.attraction)

    probs.append(outside.area)
    total_probs = sum(probs)
    probs[:] = [i / total_probs for i in probs]

    return probs, outside

# trims regional polygons so they don't include any overlap with hotspots
def trim_regpoly(reg_poly, attraction_polys):

    trimmed = []
    for reg in reg_poly:
        trim = reg[1]
        for hotspot in attraction_polys:
            trim = trim.difference(hotspot.polygon)
        trimmed.append(trim)

    trimmed_area = [trim.area for trim in trimmed]
    summed = sum(trimmed_area)
    trimmed_prob = [area/summed for area in trimmed_area]

    return trimmed_prob

# Returns array of locations give probabilities calculated from areas and attraction factors
def get_location(boundary, reg_poly, attraction_polys, fish, fish_num):

    fish_spec_polys = []

    for poly in attraction_polys:
        if poly.fish == fish:
            fish_spec_polys.append(poly)

    probs, outside = hotspot_prob(boundary, fish_spec_polys)

    outside_reg_prob = trim_regpoly(reg_poly, attraction_polys)
    in_region = np.random.choice(fish_spec_polys + [outside], size=fish_num, p=probs)

    locations = []
    for region in in_region:

        if type(region) == HotSpot:
            locations.append(region.get_chem_region())

        if type(region) == Polygon:
            locations.append(np.random.choice([row[0] for row in reg_poly], p=outside_reg_prob))


    return locations