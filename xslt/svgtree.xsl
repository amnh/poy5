<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0">
    <xsl:output method="xml"/>
    <!-- Change this parameter to false to generate branch lengths -->
    <xsl:param name="equal" select="true()"/>

    <xsl:include href="branchlength.xsl"/>

    <xsl:template name="next_length">
        <xsl:param name="node"/>
        <xsl:param name="depth"/>
        <xsl:choose>
            <xsl:when test="$equal">
                <xsl:variable name="res">
                    <xsl:call-template name="equal">
                        <xsl:with-param name="node" select="$node"/>
                        <xsl:with-param name="depth" select="$depth"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:value-of select="$res"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:variable name="res">
                    <xsl:call-template name="length">
                        <xsl:with-param name="node" select="$node"/>
                        <xsl:with-param name="depth" select="$depth"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:value-of select="$res"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

    <xsl:template match="Diagnosis">
        <forest>
        <xsl:call-template name="AllForests">
            <xsl:with-param name="list" select="/Diagnosis/Forest[*]"/>
            <xsl:with-param name="line" select="'0'"/>
        </xsl:call-template>
        </forest>
    </xsl:template>

    <xsl:template name="AllForests">
        <xsl:param name="list"/>
        <xsl:param name="line"/>
        <xsl:if test="$list">
            <xsl:variable name="newline">
                <xsl:call-template name="maxline">
                    <xsl:with-param name="node" select="$list[1]/Tree[1]"/>
                    <xsl:with-param name="line" select="$line"/>
                </xsl:call-template>
            </xsl:variable>
            <xsl:call-template name="Forest">
                <xsl:with-param name="forest" select="$list[1]"/>
                <xsl:with-param name="line" select="$line"/>
            </xsl:call-template>
            <xsl:call-template name="AllForests">
                <xsl:with-param name="list" select="$list[position () > 1]"/>
                <xsl:with-param name="line" select="$newline + 1"/>
            </xsl:call-template>
        </xsl:if>
    </xsl:template>

    <xsl:template name="Forest">
        <xsl:param name="forest"/>
        <xsl:param name="line"/>
        <xsl:variable name="ldepth">
            <xsl:call-template name="maxdepth">
                <xsl:with-param name="node" select="$forest/Tree[1]"/>
                <xsl:with-param name="depth" select="'0'"/>
            </xsl:call-template>
        </xsl:variable>
        <xsl:call-template name="prepare_tree">
            <xsl:with-param name="list" select="$forest/Tree/Tree[*]"/>
            <xsl:with-param name="name" select="$forest/Tree/Node/@Name"/>
            <xsl:with-param name="line" select="$line"/>
            <xsl:with-param name="depth" select="0"/>
            <xsl:with-param name="ldepth" select="$ldepth"/>
        </xsl:call-template>
    </xsl:template>

    <xsl:template name="maxdepth">
        <xsl:param name="node"/>
        <xsl:param name="depth"/>
        <xsl:choose>
            <xsl:when test="$node/Tree">
                <xsl:variable name="first">
                    <xsl:call-template name="maxdepth">
                        <xsl:with-param name="node" select="$node/Tree[1]"/>
                        <xsl:with-param name="depth" select="$depth + 1"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="second">
                    <xsl:call-template name="maxdepth">
                        <xsl:with-param name="node" select="$node/Tree[2]"/>
                        <xsl:with-param name="depth" select="$depth + 1"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:choose>
                    <xsl:when test="$first &gt; $second">
                        <xsl:value-of select="$first"/>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:value-of select="$second"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="$depth"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

    <xsl:template name="maxline">
        <xsl:param name="node"/>
        <xsl:param name="line"/>
        <xsl:choose>
            <xsl:when test="$node/Tree">
                <xsl:variable name="first">
                    <xsl:call-template name="maxline">
                        <xsl:with-param name="node" select="$node/Tree[1]"/>
                        <xsl:with-param name="line" select="$line"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="second">
                    <xsl:call-template name="maxline">
                        <xsl:with-param name="node" select="$node/Tree[2]"/>
                        <xsl:with-param name="line" select="$first"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:value-of select="$second"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="$line + 3"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

    <xsl:template name="prepare_tree">
        <xsl:param name="list"/>
        <xsl:param name="name"/>
        <xsl:param name="line"/>
        <xsl:param name="depth"/>
        <xsl:param name="ldepth"/>
        <xsl:choose>
            <xsl:when test="$list">
                <xsl:variable name="maxline">
                    <xsl:call-template name="maxline">
                        <xsl:with-param name="node" select="$list[2]"/>
                        <xsl:with-param name="line" select="$line"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="maxline1">
                    <xsl:call-template name="maxline">
                        <xsl:with-param name="node" select="$list[1]"/>
                        <xsl:with-param name="line" select="$maxline"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="depth_of_2">
                    <xsl:call-template name="next_length">
                        <xsl:with-param name="node" select="$list[2]/Node"/>
                        <xsl:with-param name="depth" select="$depth"/>
                    </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="depth_of_1">
                    <xsl:call-template name="next_length">
                        <xsl:with-param name="node" select="$list[1]/Node"/>
                        <xsl:with-param name="depth" select="$depth"/>
                    </xsl:call-template>
                </xsl:variable>
                <node depth="{$depth}" 
                    line="{floor (($line + $maxline1 - 1) div 2)}"
                    branchtext="{$depth}"
                    maxline="{$maxline1}">
                    <xsl:call-template name="prepare_tree">
                        <xsl:with-param name="list" select="$list[2]/Tree"/>
                        <xsl:with-param name="name" select="$list[2]/Node/@Name"/>
                        <xsl:with-param name="line" select="$line"/>
                        <xsl:with-param name="depth" select="$depth_of_2"/>
                        <xsl:with-param name="ldepth" select="$ldepth"/>
                    </xsl:call-template>
                    <xsl:call-template name="prepare_tree">
                        <xsl:with-param name="list" select="$list[1]/Tree"/>
                        <xsl:with-param name="name" select="$list[1]/Node/@Name"/>
                        <xsl:with-param name="line" select="$maxline"/>
                        <xsl:with-param name="depth" select="$depth_of_1"/>
                        <xsl:with-param name="ldepth" select="$ldepth"/>
                    </xsl:call-template>
                </node>
            </xsl:when>
            <xsl:otherwise>
                <xsl:variable name="new_ldepth">
                    <xsl:choose>
                        <xsl:when test="$equal">
                            <xsl:value-of select="$ldepth"/>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl:value-of select="$depth"/>
                        </xsl:otherwise>
                    </xsl:choose>
                </xsl:variable>
                <node depth="{$depth}" 
                    stringdepth="{$new_ldepth}"
                    line = "{$line + 2}" 
                    maxline = "{$line}" 
                    branchtext="{$depth}"
                    name="{$name}"> 
                </node>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

</xsl:transform>
