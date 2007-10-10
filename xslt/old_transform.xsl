<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="xml"/>
	
	<xsl:template match="/Forest/Tree">
	<forest>
		<tree>
			<xsl:for-each select="//Node">
				<node>	
					<!-- get node id -->
					<id>
						<xsl:value-of select="@Name"/>
					</id>
					
					<!-- transfers -->	
					
					<!-- document changes/mutations and deduce types on non-additive characters -->
					<xsl:if test="Non_Additive_Character != '' and position() != 1">
						<transformations>
						<xsl:for-each select="Non_Additive_Character">
							<xsl:variable name="characterPosition" select="position()"></xsl:variable>
							<xsl:if test="normalize-space(.) != normalize-space(../../../Node/Non_Additive_Character[$characterPosition]) and @Class='Final'">
							<!-- only for class=final; ?? potentially to be changed in the future if user needs the other sequences -->
								<transformation>
									<!-- assign variables -->
									<xsl:variable name="ancestralState" select="../../../Node/Non_Additive_Character[$characterPosition]"></xsl:variable>
									<xsl:variable name="descendantState" select="."></xsl:variable>
									<xsl:variable name="ancA" select="contains($ancestralState, '1')"></xsl:variable>
									<xsl:variable name="ancC" select="contains($ancestralState, '2')"></xsl:variable>
									<xsl:variable name="ancG" select="contains($ancestralState, '3')"></xsl:variable>
									<xsl:variable name="ancT" select="contains($ancestralState, '4')"></xsl:variable>
									<xsl:variable name="anc-" select="contains($ancestralState, '5')"></xsl:variable> <!-- gap -->
									<xsl:variable name="descA" select="contains($descendantState, '1')"></xsl:variable>
									<xsl:variable name="descC" select="contains($descendantState, '2')"></xsl:variable>
									<xsl:variable name="descG" select="contains($descendantState, '3')"></xsl:variable>
									<xsl:variable name="descT" select="contains($descendantState, '4')"></xsl:variable>
									<xsl:variable name="desc-" select="contains($descendantState, '5')"></xsl:variable>
									<xsl:variable name="ancSize">
										<xsl:choose>
											<xsl:when test="($ancA and $ancC) and ($ancG and $ancT)">
												<xsl:text>0</xsl:text>
											</xsl:when>
											<xsl:otherwise>
												<xsl:value-of select="number($ancA)+number($ancC)+number($ancG)+number($ancT)"/>
											</xsl:otherwise>
										</xsl:choose>
									</xsl:variable>
									<xsl:variable name="descSize">
										<xsl:choose>
											<xsl:when test="($descA and $descC) and ($descG and $descT)">
												<xsl:text>0</xsl:text>
											</xsl:when>
											<xsl:otherwise>
												<xsl:value-of select="number($descA)+number($descC)+number($descG)+number($descT)"/>
											</xsl:otherwise>
										</xsl:choose>
									</xsl:variable>
									<xsl:variable name="intersectionSize" select="
										number($ancA and $descA) +
										number($ancC and $descC) + 
										number($ancG and $descG) +
										number($ancT and $descT)"></xsl:variable>
										
									
									<!--<xsl:attribute name="Character">0</xsl:attribute>-->
									<xsl:attribute name="Pos">
										<xsl:value-of select="substring-after(@Name,':')"/>
									</xsl:attribute>
									<xsl:attribute name="AncS">
										<xsl:if test="$ancA">A</xsl:if>
										<xsl:if test="$ancC">C</xsl:if>
										<xsl:if test="$ancG">G</xsl:if>
										<xsl:if test="$ancT">T</xsl:if>
										<xsl:if test="$anc-">-</xsl:if>
									</xsl:attribute>
									<xsl:attribute name="DescS">
										<xsl:if test="$descA">A</xsl:if>
										<xsl:if test="$descC">C</xsl:if>
										<xsl:if test="$descG">G</xsl:if>
										<xsl:if test="$descT">T</xsl:if>
										<xsl:if test="$desc-">-</xsl:if>
									</xsl:attribute>
									
									<xsl:attribute name="Type">
										<xsl:choose>
											<xsl:when test="$ancSize = 0">
												<xsl:text>Ins</xsl:text>
											</xsl:when>
											<xsl:when test="$descSize = 0">
												<xsl:text>Del</xsl:text>
											</xsl:when>
											<xsl:otherwise>
												<xsl:choose>
													<!-- when intersection is empty set or there is a state of length 3, it's not definite -->
													<xsl:when test="$intersectionSize != 0 or ($ancSize = 3 or $descSize=3)">
														<xsl:text>ABC</xsl:text>
													</xsl:when>
													<!-- changes between states of set size 1 -->
													<xsl:when test="$ancSize = 1 and $descSize = 1">
														<xsl:choose>
															<xsl:when test="(($ancA or $ancG) and ($descC or $descT)) or
																				 (($ancC or $ancT) and ($descA or $descG))">
																<xsl:text>Tv</xsl:text>
															</xsl:when>
															<xsl:otherwise>
																<xsl:text>Ti</xsl:text>
															</xsl:otherwise>
														</xsl:choose>
													</xsl:when>
													<!-- changes between states of set size 2 -->
													<xsl:when test="$ancSize = 2 and $descSize = 2">
														<xsl:choose>
															<xsl:when test="($ancA and $ancG) or ($ancC and $ancT)">
																<xsl:text>Tv</xsl:text>
															</xsl:when>
															<xsl:otherwise>
																<xsl:text>ABC</xsl:text>
															</xsl:otherwise>
														</xsl:choose>
													</xsl:when>
													<!-- cases left: set sizes of 1 and 2 -->
													<xsl:otherwise>
														<xsl:choose>
															<xsl:when test="
																				($ancA and not($descG)) or
																				($ancC and not($descT)) or
																				($ancG and not($descA)) or
																				($ancT and not($descC))
																				">
																<xsl:text>Tv</xsl:text>
															</xsl:when>
															<xsl:otherwise>
																<xsl:text>ABC</xsl:text>
															</xsl:otherwise>
														</xsl:choose>
													</xsl:otherwise>
												</xsl:choose>
											</xsl:otherwise>
										</xsl:choose>
									</xsl:attribute>
									
									<xsl:attribute name="Cost">
										<xsl:value-of select="@Cost"/>
									</xsl:attribute>
									
									<xsl:attribute name="Definite">
										<xsl:choose>
											<xsl:when test="$intersectionSize != 0">
												<xsl:text>False</xsl:text>
											</xsl:when>
											<xsl:otherwise>
												<xsl:text>True</xsl:text>
											</xsl:otherwise>
										</xsl:choose>
									</xsl:attribute>
								</transformation>
							</xsl:if>
						</xsl:for-each>
						</transformations>
					</xsl:if>
					
					<xsl:if test="Chromosome != '' and position() != 1">
						<transformations>
						<xsl:for-each select="Chromosome">
							<xsl:variable name="characterPosition" select="position()"></xsl:variable>
							<xsl:if test="@Class='Final'">
								<transformation>
									<xsl:attribute name="AncS">
										<xsl:value-of select="ChromosomeMap/SegmentMap/@AncestorStartPosition"/>
										<xsl:text> </xsl:text>
										<xsl:value-of select="ChromosomeMap/SegmentMap/@AncestorEndPosition"/>
									</xsl:attribute>
									
									<xsl:attribute name="DescS">
										<xsl:value-of select="ChromosomeMap/SegmentMap/@DescendantStartPosition"/>
										<xsl:text> </xsl:text>
										<xsl:value-of select="ChromosomeMap/SegmentMap/@DescendantEndPosition"/>
									</xsl:attribute>
									
									<xsl:attribute name="Pos">
										<xsl:value-of select="@Reference_code"/>
									</xsl:attribute>
									
									<xsl:attribute name="Type">
										<xsl:text>???</xsl:text>
									</xsl:attribute>
									<xsl:attribute name="Cost">
										<xsl:value-of select="@Cost"/>
									</xsl:attribute>
									
									<xsl:attribute name="Definite">
										<xsl:text>???</xsl:text>
									</xsl:attribute>
								</transformation>
							</xsl:if>
						</xsl:for-each>
						</transformations>
					</xsl:if>
					
					<!-- put ancestors if not root -->
					<xsl:if test="@Name != 'root'">
						<ancestors>
							<ancestor>
									<xsl:variable name="ancestorName" select="../../Node/@Name"></xsl:variable>
								<xsl:attribute name="id">
									<xsl:value-of select="$ancestorName"/>
								</xsl:attribute>
								<xsl:attribute name="type">
									<xsl:choose>
										<xsl:when test="$ancestorName = 'root'">
											<xsl:text>root</xsl:text>
										</xsl:when>
										<xsl:otherwise>
											<xsl:text>node</xsl:text>
										</xsl:otherwise>
									</xsl:choose>
								</xsl:attribute>
							</ancestor>
						</ancestors>
					
					</xsl:if>
				</node>
			</xsl:for-each>
		</tree>
	</forest>
	</xsl:template>
	
</xsl:stylesheet>
