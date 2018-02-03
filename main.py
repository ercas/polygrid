#!/usr/bin/env python3

import fiona
import json
import math
import multiprocessing
import msgpack
import os
import shapely.geometry
import sys

FILE_FORMAT_VERSION = 0.2
THREADS = multiprocessing.cpu_count()
ROW_TEMPFILE_TEMPLATE = "polygrid_row_%d.msgpack"

def compute_row_mp(args):
    (polygrid, y) = args

    print("computing row %d (pid: %d)" % (
        y, multiprocessing.current_process().pid
    ))

    row = [[]] * polygrid.columns

    for x in range(polygrid.columns):
        row[x] = polygrid.compute_cell(x, y)

    with open(ROW_TEMPFILE_TEMPLATE % y, "wb") as f:
        msgpack.dump(row, f)

class PolyGrid(object):

    def __init__(self, rows = None, columns = None, input_file = None):
        """ Initialize new PolyGrid object

        Args:
            rows: The number of rows in the grid
            columns: The number of columns in the grid, defaults to the number
                of rows if not provided
            input_file: A file that a PolyGrid should be loaded from, if
                supplied
        """

        if (input_file is None):
            self.new_blank_grid(rows, columns)
        else:
            self.open(input_file)

    def new_blank_grid(self, rows, columns = None):
        """ Initialize a new, blank grid

        Args:
            rows: The number of rows in the grid
            columns: The number of columns in the grid, defaults to the number
                of rows if not provided
        """

        self.shapes = []

        self.bounds = None
        self.lon_range = None
        self.lat_range = None

        self.rows = rows
        if (columns == None):
            self.columns = rows
        else:
            self.columns = columns

        self.grid = [
            [
                None
                for col in range(self.columns)
            ]
            for row in range(self.rows)
        ]

    def open(self, input_file):
        """ Load a PolyGrid from a file

        Args:
            input_file: The file to load the PolyGrid from
        """

        with open(input_file, "rb") as f:
            header = f.readline()

            data = json.loads(header.decode())
            assert data["version"] == FILE_FORMAT_VERSION, "incompatible file"

            print("loading shapes")
            self.shapes = [
                {
                    "geoid": shape["geoid"],
                    "geometry": shapely.geometry.shape(shape["geometry"])
                }
                for shape in data["shapes"]
            ]
            self.bounds = data["bounds"]
            self.lon_range = self.bounds[2] - self.bounds[0]
            self.lat_range = self.bounds[3] - self.bounds[1]

            """
            self.grid = []
            for row in f:
                row_list = []
                for cell in row.split(" "):
                    try:
                        row_list.append(msgpack.dumps(
                            [
                            int(i)
                            for i in cell.split(",")
                            ]
                        ))
                    except ValueError:
                        row_list.append([])
                self.grid.append(row_list)
            """

            print("loading grid")
            self.grid = msgpack.load(f)
            self.columns = len(self.grid[0])
            self.rows = len(self.grid)

    def save(self, output_file):
        """ Save a PolyGrid to a file

        Args:
            output_file: The file to save the PolyGrid to
        """

        with open(output_file, "wb") as f:
            f.write(
                bytes(
                    "%s\n" % json.dumps(
                        {
                            "version": FILE_FORMAT_VERSION,
                            "shapes": [
                                {
                                    "geoid": shape["geoid"],
                                    "geometry": shapely.geometry.mapping(
                                        shape["geometry"]
                                    )
                                }
                                for shape in self.shapes
                            ],
                            "bounds": self.bounds
                        }
                    ),
                    "utf-8"
                )
            )
            msgpack.dump(self.grid, f)
            """
            for row in self.grid:
                f.write("%s\n" %
                    " ".join([
                        ",".join([
                            str(i)
                            for i in cell
                        ])
                        for cell in row
                    ])
                )
            """

        print("wrote to %s" % output_file)

    def to_grid(self, lon, lat):
        """ Convert a longitude, latitude pair to a coordinate pair that points
        to a grid cell

        Args:
            lon: The longitude to be used
            lat: The latitude to be used

        Returns:
            A tuple containing the x and y coordinates of the grid cell
        """


        return (
            math.floor((lon - self.bounds[0]) / self.lon_range * self.columns),
            math.floor((lat - self.bounds[1]) / self.lat_range * self.rows)
        )

    def to_world(self, x, y):
        """ Convert a coordinate pair that points to a grid cell to a
        longitude, latitude pair

        Args:
            x: The column of the grid cell
            y: The row of the grid cell

        Returns:
            A tuple containing the longitude and latitude equivalents of the
            grid cell's lower left corner
        """

        return (
            self.bounds[0] + (x / self.columns * self.lon_range),
            self.bounds[1] + (y / self.rows * self.lat_range)
        )

    def compute_cell(self, x, y):
        """ Fill in a single cell of the grid

        Args:
            x: The column of the grid cell
            y: The row of the grid cell

        Returns:
            The contents of the cell at that position
        """

        southwest = self.to_world(x, y)
        northeast = self.to_world(x + 1, y + 1)

        cell_poly = shapely.geometry.box(
            southwest[0], southwest[1], northeast[0], northeast[1]
        )

        cell = []
        for shape_index in range(len(self.shapes)):
            if (self.shapes[shape_index]["geometry"].intersects(cell_poly)):
                cell.append(shape_index)

        #return cell
        return msgpack.dumps(cell)

    def compute_row(self, y):
        """ Fill in a single row of the grid

        Args:
            y: The row of the grid cell
        """

        print("computing row %d (pid: %d)" % (
            y, multiprocessing.current_process().pid
        ))
        for x in range(self.columns):
            self.grid[y][x] = self.compute_cell(x, y)

    def recompute(self):
        """ Recompute the extents of all shapes and the polygrid

        This function must be run when shapes are added so that lookup
        functions reflect the current state of the PolyGrid.
        """

        ## bounds recomputation ################################################
        bboxes = [shape["geometry"].bounds for shape in self.shapes]

        self.bounds = [
            min([bbox[i] for bbox in bboxes])
                if i <= 1
                else max([bbox[i] for bbox in bboxes])
            for i in range(4)
        ]
        self.lon_range = self.bounds[2] - self.bounds[0]
        self.lat_range = self.bounds[3] - self.bounds[1]

        print("recomputed bounds for %d shapes" % len(self.shapes))

        ## grid cell intersection recomputation ################################
        self.grid = [[]] * self.rows
        #self.grid = [[ [[]] * self.columns ]] * self.rows

        pool = multiprocessing.Pool(THREADS)
        pool.map(
            compute_row_mp,
            zip(
                [self] * self.rows,
                range(self.rows)
            )
        )
        pool.close()
        pool.join()

        print("merging rows")
        for y in range(self.rows):
            row_path = ROW_TEMPFILE_TEMPLATE % y
            with open(row_path, "rb") as f:
                self.grid[y] = msgpack.load(f)
            os.remove(row_path)

    def add_shapes(self, filename):
        """ Add shapes to the grid

        Args:
            filename: The path to a shapefile, GeoJSON, or other file
                containing shapes supported by the fiona library
        """

        i = 0
        with fiona.drivers(), fiona.open(filename) as src:
            for shape in src:
                shape["geometry"] = shapely.geometry.shape(shape.pop("geometry"))
                shape["geoid"] = shape.pop("properties")["GEOID"]
                del shape["id"]
                del shape["type"]
                self.shapes.append(shape)

                i += 1
                sys.stdout.write("\rloaded shape %d" % i)
                sys.stdout.flush()
            print("added %d shapes" % len(src))

    def lookup(self, lon, lat):
        grid_coords = self.to_grid(lon, lat)

        try:
            #cell = self.grid[grid_coords[1]][grid_coords[0]]
            cell = msgpack.loads(self.grid[grid_coords[1]][grid_coords[0]])
            if (len(cell) == 1):
                return self.shapes[cell[0]]["geoid"]
            elif (len(cell) == 0):
                return None
            else:
                point = shapely.geometry.Point(lon, lat)

                for shape_index in cell:
                    shape = self.shapes[shape_index]
                    if (shape["geometry"].contains(point)):
                        return shape["geoid"]

        except IndexError:
            pass
