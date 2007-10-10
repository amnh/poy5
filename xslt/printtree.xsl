<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0">
    <xsl:output method="xml" doctype-public="-//W3C//DTD SVG 1.1//EN"
        doctype-system="http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd"/>

    <!-- A transform to produce a tree in svg format. -->

    <!-- USER PASSED (OR DEFINED) PARAMETERS -->
    <!-- Overall spacing between the vertices of the tree, in the vertical and
    horizontal coordinates -->
    <xsl:param name="y_space" select="'5'"/>
    <xsl:param name="x_space" select="'10'"/>
    <!-- The stroke width -->
    <xsl:param name="width" select="'1'"/>
    <xsl:param name="print_branchtext" select="false()"/>
    <xsl:param name="hl_terminal1" select="'first'"/>
    <xsl:param name="hl_terminal2" select="'second'"/>

    <!-- END OF USER PARAMENTERS -->

    <!-- A rule to match a forest of trees to be presented. -->
    <xsl:template match="forest">
        <svg width="100%" height="100%" version="1.1">
            <xsl:for-each select="node">
                <xsl:call-template name="print_tree">
                    <xsl:with-param name="tree" select="."/>
                </xsl:call-template>
            </xsl:for-each>
        </svg>
    </xsl:template>

    <!-- A function to count how many of the desired taxa are in this clade -->

    <!-- The recursive tree printer in svg. -->
    <xsl:template name="print_tree">
        <xsl:param name="tree"/>
        <xsl:variable name="color">
            <xsl:choose>
                <xsl:when test="$tree/@color">stroke:'<xsl:value-of select="$tree/@color"/>'</xsl:when>
                <xsl:otherwise>stroke:rgb(0,0,0)</xsl:otherwise>
        </xsl:choose>
        </xsl:variable>
        <xsl:choose>
            <xsl:when test="$tree/*">
                <!-- We print the vertical line -->
                <line x1="{$tree/@depth * $x_space}" 
                    y1="{$tree/node[1]/@line * $y_space}" 
                    x2="{$tree/@depth * $x_space}" 
                    y2="{$tree/node[2]/@line * $y_space}"
                    style="{$color}; stroke-width:{$width}"/>
                <!-- We print the upper horizontal line of the subtree -->
                <line x1="{$tree/@depth * $x_space}" 
                    y1="{$tree/node[1]/@line * $y_space}"
                    x2="{$tree/node[1]/@depth * $x_space}" 
                    y2 = "{$tree/node[1]/@line * $y_space}"
                    style="{$color}; stroke-width:{$width}"/>
                <!-- We print the lower horizontal line of the subtree -->
                <line x1="{$tree/@depth * $x_space}" 
                    y1="{$tree/node[2]/@line * $y_space}"
                    x2="{$tree/node[2]/@depth * $x_space}" 
                    y2 = "{$tree/node[2]/@line * $y_space}"
                    style="{$color}; stroke-width:{$width}"/>
                <!-- We decide if we print or not the the branches lengths below -->
                <xsl:if test="$print_branchtext">
                    <text x="{($tree/@depth + (floor (($tree/node[2]/@depth - $tree/@depth) div 2))) * $x_space}" y="{($tree/node[2]/@line * $y_space) - 1}">
                        <xsl:value-of select="$tree/node[2]/@branchtext - $tree/@branchtext"/>
                    </text>
                    <text x="{($tree/@depth + (floor(($tree/node[1]/@depth - $tree/@depth) div 2))) * $x_space}" y="{($tree/node[1]/@line * $y_space) - 1}">
                        <xsl:value-of select="$tree/node[1]/@branchtext - $tree/@branchtext"/>
                    </text>
                </xsl:if>

                <xsl:call-template name="print_tree">
                    <xsl:with-param name="tree" select="$tree/node[1]"/>
                </xsl:call-template>
                <xsl:call-template name="print_tree">
                    <xsl:with-param name="tree" select="$tree/node[2]"/>
                </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
                <!-- We print the string of the leaf -->
                <line x1="{$tree/@depth * $x_space}" y1 = "{$tree/@line * $y_space}"
                    x2="{$tree/@stringdepth * $x_space}" y2 = "{$tree/@line * $y_space}"/>
                <text x="{$tree/@stringdepth * $x_space}" y="{$tree/@line * $y_space}">
                    <xsl:value-of select="$tree/@name"/>
                </text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

</xsl:transform>
