<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0" xmlns:phyc="http://commandLanguage.phycas.org">
<xsl:output method="xml" indent="yes"/>
<xsl:strip-space elements="phyc:command_list"/>
<xsl:template match="/">

<map version="1.0">
	<mapID target="intro" url="master/intro.html" />
	<mapID target="usingGUI" url="master/intro.html#GUI" />
	<mapID target="commands" url="master/commands.html" />
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:sort select="@gui_label" data-type="number"/>
                            <xsl:call-template name="command_map"/>
                        </xsl:for-each>
                    </xsl:when>
				</xsl:choose>
			</xsl:for-each>
</map>

</xsl:template>

<xsl:template name="command_map">
	<xsl:variable name="command_label"><xsl:value-of select="./@label"/></xsl:variable>
    <mapID url="master/commands.html#command{$command_label}" target="command{$command_label}"/>
		<xsl:for-each select="cmd_param_group">
			<xsl:for-each select="cmd_param">
				<mapID url="master/commands.html#command{$command_label}Param{./@label}{./@placement}" target="command{$command_label}Param{./@label}{./@placement}"/>
			</xsl:for-each>
		</xsl:for-each>
</xsl:template>


</xsl:stylesheet>
