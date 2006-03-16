<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0" xmlns:phyc="http://commandLanguage.phycas.org">
<xsl:output method="html" indent="yes"/>
<xsl:strip-space elements="phyc:command_list"/>
<xsl:template match="/">

<HTML>
<HEAD>
<TITLE>Phycas Commands</TITLE>
</HEAD>
<BODY BGCOLOR="#ffffff">
<H1>Phycas Commands</H1>
<span style="font-size: 12.0pt; font-family: Times New Roman">The Phycas Software consists of a series of commands that can be executed.
Each command and it's options are explained below.</span><p>
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:sort select="@gui_label" data-type="number"/>
                            <xsl:call-template name="command_index"/>
                        </xsl:for-each>
                    </xsl:when>
				</xsl:choose>
			</xsl:for-each>

</p>

</BODY></HTML>

</xsl:template>

<xsl:template name="command_index">
	<xsl:variable name="command_label"><xsl:value-of select="./@label"/></xsl:variable>
	<xsl:variable name="command_description"><xsl:value-of select="description"/></xsl:variable>
	<font size="5"><a name="command{$command_label}"><xsl:value-of select="$command_label"/></a></font>
	<p><xsl:value-of select="$command_description"/></p>
	<p><b>Command Options:</b>
	<ul>
		<xsl:for-each select="cmd_param_group">
			<xsl:for-each select="cmd_param">
				<li>
				<font size="4"><a name="command{$command_label}Param{./@label}{./@placement}"><xsl:value-of select="./gui_description"/></a></font>
				<p><xsl:value-of select="./description"/></p> <!--Probably need another field for help info-->
				<p><i>GUI:  </i>  
				<xsl:choose>
					<xsl:when test="./mixed_type_info">
							<xsl:for-each select="./mixed_type_info/cmd_param">
									<xsl:choose>
										<xsl:when test="position()=last()">
											<xsl:call-template name="process_type_info">
													<xsl:with-param name="mixed">false</xsl:with-param>
											</xsl:call-template>
										</xsl:when>
										<xsl:otherwise>
											<xsl:call-template name="process_type_info">
													<xsl:with-param name="mixed">true</xsl:with-param>
											</xsl:call-template>
										</xsl:otherwise>
									</xsl:choose>
							</xsl:for-each>	
					</xsl:when>
					<xsl:otherwise>
						<xsl:call-template name="process_type_info"/>
					</xsl:otherwise>
				</xsl:choose>
				</p></li>
			</xsl:for-each>
		</xsl:for-each>
	</ul></p>	
</xsl:template>

<xsl:template name="process_type_info">
	<xsl:param name="mixed">false</xsl:param>
				<xsl:if test="./bool_type_info">
					Select the checkbox for this option to be true, and deselect for false
				</xsl:if>
				<xsl:if test="./char_set_type_info or ./tax_set_type_info or ./tree_set_type_info">
					To use a set that already exists, select it from the drop-down.  Items can then be added to or removed from the set by selecting or deselecting them from the list above the drop-down.
					The text field below the drop-down can be used to directly type in a set of items, instead of having to select them from the list.
				</xsl:if>
				<xsl:if test="./choice_type_info">
					Either select a radio button or an item from the list
				</xsl:if>
				<xsl:if test="./integer_type_info">
					Enter an integer value into the text field
				</xsl:if>
				<xsl:if test="./double_type_info">
					Enter a number (whole or decimal) into the text field
				</xsl:if>
				<xsl:if test="./infile_type_info">
					Enter the name of an input file into the text field, or select the choose button to browse to the file location
				</xsl:if>
				<xsl:if test="./name_type_info or ./restricted_string_type_info">
					Enter a value into the text field.  This option has some values that are not allowed.  If one of these is entered, an error message will be displayed in the message box
				</xsl:if>
				<xsl:if test="./output_type_info">
					Select the file checkbox to have output saved to a file.  Then enter the file name of the output file or select the choose button to browse to the file, and select a radio button to 
					either append to the file or replace it.  Select the output checkbox to have the output displayed in the Output text of the Output Window.  If available, select the plot checkbox to
					create a traceplot of the particular values.
				</xsl:if>
				<xsl:if test="./string_type_info">
					Enter a value into the text field
				</xsl:if>
				<xsl:if test="./distribution_type_info">
					Enter a probability distribution into the text field or select the choose button for help in creating one.
				</xsl:if>
				<xsl:if test="$mixed = 'true'">OR</xsl:if>
				<xsl:if test="$mixed = 'false'">.</xsl:if>
</xsl:template>
</xsl:stylesheet>
