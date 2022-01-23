/*
 * Copyright (c) 2019 Gabriel Roldan
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
package org.locationtech.jts.io.twkb;

import java.util.Objects;
import java.util.function.Function;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

/**
 * <pre>
 * {@code
 * 
 * twkb                    := <header> <geometry_body>
 * header                  := <type_and_precision> <metadata_header> [extended_dimensions_header] [geometry_body_size]
 * type_and_precision      := byte := <type_mask OR precision>)
 * type_mask               := <ubyte> (0b0000XXXX -> 1=point, 2=linestring, 3=polygon, 4=multipoint, 
 *                                     5=multilinestring, 6=multipolygon, 7=geometry collection)
 * precision               := <signed byte> (zig-zag encoded 4-bit signed integer, 0bXXXX0000. Number of base-10 decimal places 
 *                                           stored. A positive retaining information to the right of the decimal place, negative 
 *                                           rounding up to the left of the decimal place)  
 * metadata_header := byte := <bbox_flag OR size_flag OR idlist_flag OR extended_precision_flag OR empty_geometry_flag>
 * bbox_flag               := 0b00000001
 * size_flag               := 0b00000010
 * idlist_flag             := 0b00000100
 * extended_precision_flag := 0b00001000
 * empty_geometry_flag     := 0b00010000
 * 
 * # extended_dimensions_header present iif extended_precision_flag == 1
 * extended_dimensions_header  := byte := <Z_dimension_presence_flag OR M_dimension_presence_flag OR Z_precision OR M_precision>
 * Z_dimension_presence_flag   := 0b00000001 
 * M_dimension_presence_flag   := 0b00000010
 * Z_precision                 := 0b000XXX00 3-bit unsigned integer using bits 3-5 
 * M_precision                 := 0bXXX00000 3-bit unsigned integer using bits 6-8
 * 
 * # geometry_body_size present iif size_flag == 1 
 * geometry_body_size := uint32 # size in bytes of <geometry_body>
 * 
 * # geometry_body present iif empty_geometry_flag == 0
 * geometry_body := [bounds] [idlist] <geometry>
 * # bounds present iff bbox_flag == 1 
 * # 2 signed varints per dimension. i.e.:
 * # [xmin, deltax, ymin, deltay]                              iif Z_dimension_presence_flag == 0 AND M_dimension_presence_flag == 0
 * # [xmin, deltax, ymin, deltay, zmin, deltaz]                iif Z_dimension_presence_flag == 1 AND M_dimension_presence_flag == 0
 * # [xmin, deltax, ymin, deltay, zmin, deltaz, mmin, deltam]  iif Z_dimension_presence_flag == 1 AND M_dimension_presence_flag == 1
 * # [xmin, deltax, ymin, deltay, mmin, deltam]                iif Z_dimension_presence_flag == 0 AND M_dimension_presence_flag == 1
 * bounds          := sint32[4] | sint32[6] | sint32[8] 
 * geometry        := point | linestring | polygon | multipoint | multilinestring | multipolygon | geomcollection
 * point           := sint32[dimension]
 * linestring      := <npoints:uint32> [point[npoints]]
 * polygon         := <nrings:uint32> [linestring]
 * multipoint      := <nmembers:uint32> [idlist:<sint32[nmembers]>] [point[nmembers]]
 * multilinestring := <nmembers:uint32> [idlist:<sint32[nmembers]>] [linestring[nmembers]]
 * multipolygon    := <nmembers:uint32> [idlist:<sint32[nmembers]>] [polygon[nmembers]]
 * geomcollection  := <nmembers:uint32> [idlist:<sint32[nmembers]>] [twkb[nmembers]]
 * 
 * uint32 := <Unsigned variable-length encoded integer>
 * sint32 := <Signed variable-length, zig-zag encoded integer>
 * byte := <Single octect>
 * 
 * }
 * </pre>
 */
class TWKBHeader {

    public TWKBHeader() {

    }

    public TWKBHeader(TWKBHeader other) {
        this.geometryType = other.geometryType;
        this.xyPrecision = other.xyPrecision;
        this.hasBBOX = other.hasBBOX;
        this.hasSize = other.hasSize;
        this.hasIdList = other.hasIdList;
        this.isEmpty = other.isEmpty;
        this.hasZ = other.hasZ;
        this.hasM = other.hasM;
        this.zPrecision = other.zPrecision;
        this.mPrecision = other.mPrecision;
        this.geometryBodySize = other.geometryBodySize;
    }

    public GeometryType geometryType() {
        return this.geometryType;
    }

    public int xyPrecision() {
        return this.xyPrecision;
    }

    public boolean hasBBOX() {
        return this.hasBBOX;
    }

    public boolean hasSize() {
        return this.hasSize;
    }

    public boolean hasIdList() {
        return this.hasIdList;
    }

    public boolean isEmpty() {
        return this.isEmpty;
    }

    public boolean hasZ() {
        return this.hasZ;
    }

    public boolean hasM() {
        return this.hasM;
    }

    public int zPrecision() {
        return this.zPrecision;
    }

    public int mPrecision() {
        return this.mPrecision;
    }

    public TWKBHeader setGeometryType(GeometryType geometryType) {
        this.geometryType = geometryType;
        return this;
    }

    public TWKBHeader setXyPrecision(int xyPrecision) {
        this.xyPrecision = xyPrecision;
        return this;
    }

    public TWKBHeader setHasBBOX(boolean hasBBOX) {
        this.hasBBOX = hasBBOX;
        return this;
    }

    public TWKBHeader setHasSize(boolean hasSize) {
        this.hasSize = hasSize;
        return this;
    }

    public TWKBHeader setHasIdList(boolean hasIdList) {
        this.hasIdList = hasIdList;
        return this;
    }

    public TWKBHeader setEmpty(boolean empty) {
        isEmpty = empty;
        return this;
    }

    public TWKBHeader setHasZ(boolean hasZ) {
        this.hasZ = hasZ;
        return this;
    }

    public TWKBHeader setHasM(boolean hasM) {
        this.hasM = hasM;
        return this;
    }

    public TWKBHeader setZPrecision(int zPrecision) {
        this.zPrecision = zPrecision;
        return this;
    }

    public TWKBHeader setMPrecision(int mPrecision) {
        this.mPrecision = mPrecision;
        return this;
    }

    public TWKBHeader setGeometryBodySize(int geometryBodySize) {
        this.geometryBodySize = geometryBodySize;
        return this;
    }

    @Override
    public String toString() {
        return "TWKBHeader{" +
            "geometryType=" + geometryType +
            ", xyPrecision=" + xyPrecision +
            ", hasBBOX=" + hasBBOX +
            ", hasSize=" + hasSize +
            ", hasIdList=" + hasIdList +
            ", isEmpty=" + isEmpty +
            ", hasZ=" + hasZ +
            ", hasM=" + hasM +
            ", zPrecision=" + zPrecision +
            ", mPrecision=" + mPrecision +
            ", geometryBodySize=" + geometryBodySize +
            '}';
    }

    public int geometryBodySize() {
        return this.geometryBodySize;
    }

    enum GeometryType {
        POINT(1, GeometryFactory::createPoint), //
        LINESTRING(2, GeometryFactory::createLineString), //
        POLYGON(3, GeometryFactory::createPolygon), //
        MULTIPOINT(4, GeometryFactory::createMultiPoint), //
        MULTILINESTRING(5, GeometryFactory::createMultiLineString), //
        MULTIPOLYGON(6, GeometryFactory::createMultiPolygon), //
        GEOMETRYCOLLECTION(7, GeometryFactory::createGeometryCollection);

        private final int value;

        private final Function<GeometryFactory, Geometry> emptyBuilder;

        GeometryType(int value, Function<GeometryFactory, Geometry> emptyBuilder) {
            this.value = value;
            this.emptyBuilder = emptyBuilder;
        }

        public int getValue() {
            return value;
        }

        // held as a class variable cause calling values() on a tight loop has a non-depreciable
        // performance impact
        private static final GeometryType[] VALUES = GeometryType.values();

        public static GeometryType valueOf(int value) {
            return VALUES[value - 1];
        }

        public static GeometryType valueOf(Class<? extends Geometry> gclass) {
            Objects.requireNonNull(gclass);
            if (Point.class.isAssignableFrom(gclass))
                return POINT;
            if (LineString.class.isAssignableFrom(gclass))
                return LINESTRING;
            if (Polygon.class.isAssignableFrom(gclass))
                return POLYGON;
            if (MultiPoint.class.isAssignableFrom(gclass))
                return MULTIPOINT;
            if (MultiLineString.class.isAssignableFrom(gclass))
                return MULTILINESTRING;
            if (MultiPolygon.class.isAssignableFrom(gclass))
                return MULTIPOLYGON;
            if (GeometryCollection.class.isAssignableFrom(gclass))
                return GEOMETRYCOLLECTION;

            throw new IllegalArgumentException("Unrecognized geometry tpye: " + gclass);
        }

        public Geometry createEmpty(GeometryFactory factory) {
            return this.emptyBuilder.apply(factory);
        }
    }

    // first 1-byte header //
    private GeometryType geometryType;

    private int xyPrecision = 0;

    // metadata_header := byte
    // bbox_flag := 0b00000001
    // size_flag := 0b00000010
    // idlist_flag := 0b00000100
    // extended_precision_flag := 0b00001000
    // empty_geometry_flag := 0b00010000
    private boolean hasBBOX = false;

    private boolean hasSize = false;

    private boolean hasIdList = false;

    public boolean hasExtendedPrecision() {
        return hasZ() || hasM();
    }

    private boolean isEmpty = false;

    // extended_dimensions_header present iif extended_precision_flag == 1
    // extended_dimensions_header := byte
    // Z_dimension_presence_flag := 0b00000001
    // M_dimension_presence_flag := 0b00000010
    // Z_precision := 0b000XXX00 3-bit unsigned integer using bits 3-5
    // M_precision := 0bXXX00000 3-bit unsigned integer using bits 6-8

    private boolean hasZ = false;

    private boolean hasM = false;

    private int zPrecision = 0;

    private int mPrecision = 0;

    /**
     * Size of encoded geometry body, iif size_flag == 1, defaults to {@code -1} if {@link #hasSize}
     * ({@code == false}
     * <p>
     * {@code geometry_body_size := uint32 # size in bytes of <geometry_body>}
     */
    private int geometryBodySize;

    /////////////////////// custom optimizations ///////////////////////

    public int getDimensions() {
        return 2 + (hasZ ? 1 : 0) + (hasM ? 1 : 0);
    }

    public int getPrecision(final int dimensionIndex) {
        switch (dimensionIndex) {
        case 0:
        case 1:
            return xyPrecision;
        case 2:
            if (!(hasZ || hasM)) {
                throw new IllegalArgumentException("Geometry only has XY dimensions.");
            }
            return hasZ ? zPrecision : mPrecision;
        case 3:
            if (!(hasZ && hasM)) {
                throw new IllegalArgumentException("Geometry has no M dimension.");
            }
            return mPrecision;
        default:
            throw new IllegalArgumentException(
                    "Dimension index shall be between 0 and 3: " + dimensionIndex);
        }
    }

}