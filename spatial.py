from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from copy import deepcopy

class HotSpot:

    # Potential point of confusion:
    # When hot spots are defined by fraction, these are fractions of regional areas
    # NOT  what fraction of the hotspot is region x.

    def __init__(self, name, deftype, assos_fish, attraction, list):

        self.name = name
        self.deftype = deftype
        self.fish = assos_fish
        self.attraction = attraction
        self.definition = list
        self.parse_list()
        self.polygon = None
        self.weights = []
        self.area = 0

    def parse_list(self):

        if self.deftype == 'Polygon':

            for i in range(len(self.definition)):
                self.definition[i] = [float(j) for j in self.definition[i].replace(' ', '').split(',')]


    def addpoly(self, polygon):

        self.polygon = polygon

    def calcweights(self, reg_polys):
        if self.deftype == 'Polygon':
            for poly in reg_polys:
                area_of_intersect = poly[1].intersection(self.polygon).area
                self.weights.append([poly[0], area_of_intersect])
        else:
            for i in range(len(reg_polys)):
                area_reg = reg_polys[i][1].area
                frac_reg = self.definition[i]
                part_area = area_reg*frac_reg
                self.area += part_area
                self.weights.append([reg_polys[i][0], part_area])

        weight_sum = sum(row[1] for row in self.weights)

        for i in range(len(self.weights)):
            unnormed_prob = self.weights[i][1]
            self.weights[i][1] = unnormed_prob / weight_sum

    def get_chem_region(self):

        names = [row[0] for row in self.weights]
        prob = [row[1] for row in self.weights]

        return np.random.choice([i for i in range(len(prob))], p=prob)


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
        #print(hotspot.polygon)
        if hotspot.polygon is not None:
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
            polygon = Polygon(hotspot.definition)
            hotspot.addpoly(polygon)
            hotspot.calcweights(reg_polygons)
        else:
            hotspot.calcweights(reg_polygons)


    # Use to see regions
    #plot_vor(vor, boundary, attract_poly)

    return boundary, reg_polygons, attract_poly

# returns probability of being in each hotspot, and probability of being outside of hotspots
def hotspot_prob(boundary, attraction_polys, regions):

    probs = []

    if attraction_polys[0].deftype == 'Polygon':
        outside = boundary
        for poly in attraction_polys:
            outside = outside.difference(poly.polygon)
            attraction_factor = poly.attraction
            probs.append(poly.polygon.area * attraction_factor)

        probs.append(outside.area)

    else:
        outside = 'Out'
        bound_area = boundary.area
        for poly in attraction_polys:
            area = poly.area
            probs.append(area * poly.attraction)

        prob_outside = bound_area - sum(probs)
        if prob_outside > 0:
            probs.append(bound_area - sum(probs))
        else:
            probs.append(0)

    total_probs = sum(probs)
    probs[:] = [i / total_probs for i in probs]


    return probs, outside

# trims regional polygons so they don't include any overlap with hotspots
def trim_regpoly(reg_poly, attraction_polys):

    trimmed = []

    if attraction_polys[0].deftype == 'Polygon':

        for reg in reg_poly:
            trim = reg[1]
            for hotspot in attraction_polys:
                trim = trim.difference(hotspot.polygon)
            trimmed.append(trim)

        trimmed_area = [trim.area for trim in trimmed]

        summed = sum(trimmed_area)
        trimmed_prob = [area / summed for area in trimmed_area]

    else:

        trimmed_prob = []

        for i in range(len(reg_poly)):
            cur_area = reg_poly[i][1].area
            total_area = reg_poly[i][1].area
            for hotspot in attraction_polys:
                frac_m = hotspot.definition[i]
                cur_area -= (frac_m * total_area)
            trimmed_prob.append(cur_area)

        summed = sum(trimmed_prob)
        trimmed_normed_prob = [area / summed for area in trimmed_prob]
        trimmed_prob = trimmed_normed_prob


    return trimmed_prob

# Returns array of locations give probabilities calculated from areas and attraction factors
def location_step(boundary, reg_poly, attraction_polys, fish, draw_num):

    fish_spec_polys = []

    for poly in attraction_polys:
        if poly.fish == fish:
            fish_spec_polys.append(poly)
    probs, outside = hotspot_prob(boundary, fish_spec_polys, reg_poly)

    outside_reg_prob = trim_regpoly(reg_poly, attraction_polys)


    return fish_spec_polys + [outside], probs, [ i for i in range(len(reg_poly))], outside_reg_prob


def new_draw(hotspotnames, probs, regnames, outside_reg_prob, draw_num):

    in_region = np.random.choice(hotspotnames, size=draw_num, p=probs)

    locations = []
    for region in in_region:

        if type(region) == HotSpot:
            locations.append(region.get_chem_region())

        if type(region) == Polygon or region == 'Out':
            locations.append(np.random.choice(regnames, p=outside_reg_prob))

    return locations


def adjust_diet_to_region(fishname, fish_loc, diet_data, fish_by_region, number_others):

    new_diet = deepcopy(diet_data[fishname])

    fish_in_region = fish_by_region[fish_loc]
    for i in range(number_others + 1, len(new_diet)):
        if new_diet[i][0] not in fish_in_region:
            new_diet[i][1] = 0


    fracs = [row[1] for row in new_diet]
    sum_frac = sum(fracs)

    if sum_frac != 0:
        for entry in new_diet:
            entry[1] = entry[1]/sum_frac

    return new_diet