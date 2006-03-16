<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:phyc="http://commandLanguage.phycas.org" version="1.0">
  <xsl:output method="xml" indent="yes"/>
<xsl:strip-space elements="command_list"/>
<xsl:template match="/">
<doc>
         <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:apply-templates select="."/>
                        </xsl:for-each>
                    </xsl:when>
                    <xsl:otherwise><xsl:for-each select="./command">
                            <xsl:if test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                                 <xsl:apply-templates select="."/>
                            </xsl:if>
                        </xsl:for-each></xsl:otherwise>
                </xsl:choose>
             </xsl:for-each>
       
</doc>
</xsl:template>

<xsl:template match="command">
    <panel id="PNL_{./@label}" focusable="false">
        <vbox alignmentX="0.5"  focusable="false">
			<label id="lbl_{./@label}" Font="Serif-Bold-18" Foreground="black" text="{./@label}" alignmentX="0.5" focusable="false"/>
			<label id="lbl_{./@label}_{@gui_order}_vboxSpace1" text="	" focusable="false"/>
			<textarea id="taDesc_{./@label}" columns = "50" Font="Serif-14" editable="false" linewrap="true" wrapstyleword ="true" Foreground="black" background="grey" text="{description}" alignmentX="0.5" focusable="false"/>
			<label id="lbl_{./@label}_{@gui_order}_vboxSpace2" text="	" focusable="false"/>
			<xsl:if test = "string(./cmd_param_group)"> 
				<xsl:for-each select="cmd_param_group">
					<xsl:sort select="@gui_order" data-type="number"/>
					<vbox border="TitledBorder({@name})" alignmentX="0.5" focusable="false">
					<xsl:for-each select="cmd_param">
						<xsl:sort select="string-length(./gui_description)" data-type="number"/>
						<xsl:if test="position()=last()">
							<xsl:variable name="gui_desc_length"><xsl:value-of select="string-length(./gui_description)"/></xsl:variable>
							<!--xsl:text>GUI Length = </xsl:text><xsl:value-of select="$gui_desc_length"></xsl:value-of-->
							<xsl:for-each select="../cmd_param">
								<xsl:sort select="@gui_order" data-type="number"/>
								<xsl:call-template name="process_cmd_param">
									<xsl:with-param name="cmd_name"><xsl:value-of select="../../@label"/></xsl:with-param>
									<xsl:with-param name="max_gui_desc_length"><xsl:value-of select="$gui_desc_length"/></xsl:with-param>
								</xsl:call-template>
							</xsl:for-each>
						</xsl:if>
					</xsl:for-each>
					</vbox>
					<hbox alignmentX="0.5" focusable="false">
						<label id="lbl_{../@label}_{@gui_order}_vboxSpace3" text="	" focusable="false"/>
					</hbox>
				</xsl:for-each>
			</xsl:if>
			<hbox alignmentX="0.5" focusable="false">
				<button text="OK" Action="okAction_{./@label}" ActionCommand = "okAction"/>
			</hbox>
        </vbox>
    </panel>
</xsl:template>

<xsl:template name="process_cmd_param">
	<xsl:param name="mixed">false</xsl:param>
	<xsl:param name="cmd_opt_desc"><xsl:value-of select="description"/></xsl:param>
	<xsl:param name="cmd_opt_gui_desc"><xsl:value-of select="gui_description"/></xsl:param>
	<xsl:param name="max_gui_desc_length">0</xsl:param>
	<xsl:param name="cmd_name">default</xsl:param>
	<xsl:variable name="cmd_name_"><xsl:call-template name="add_cmd_name"><xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param></xsl:call-template></xsl:variable>
    <xsl:variable name="cmd_opt_label">lbl_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_label_out">lbl_Out_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_label_space">lbl_<xsl:value-of select="$cmd_name_"/>_space</xsl:variable>
    <xsl:variable name="cmd_opt_text_field_string">tf_String_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_text_field_res_string">tf_RestrictedString_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_text_field_file">tf_File_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_text_field_set">tf_Set_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_text_field_out">tf_Out_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_text_field_number">tf_Number_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_checkbox">cb_Bool_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_checkbox_out">cb_Out_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_combobox_class">phycasGUI.swixml.model.ComboModel_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_combobox">cmb_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_combobox_choice">cmb_Choice_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_combobox_set">cmb_Set_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_combobox_action">action_cmb_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_list_class">phycasGUI.swixml.model.ListModel_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_list">list_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_list_res_string">list_RestrictedString_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_list_set">list_Set_<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_list_label"><xsl:value-of select="gui_description"/> cannot be one of these:</xsl:variable>
    <xsl:variable name="cmd_opt_list_action">action_list<xsl:value-of select="$cmd_name_"/></xsl:variable>
    <!--xsl:variable name="cmd_opt_gui_desc"><xsl:value-of select="gui_description"/></xsl:variable-->
    <xsl:variable name="cmd_opt_choose"><xsl:value-of select="$cmd_name_"/></xsl:variable>
    <xsl:variable name="cmd_opt_default"><xsl:choose>
        <xsl:when test="bool_type_info"><xsl:value-of  select="bool_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="char_set_type_info"><xsl:value-of select="char_set_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="choice_type_info"><xsl:value-of select="choice_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="double_type_info"><xsl:value-of select="double_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="infile_type_info"><xsl:value-of select="infile_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="outfile_type_info"><xsl:value-of select="outfile_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="restricted_string_type_info"><xsl:value-of select="restricted_string_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="name_type_info"><xsl:value-of select="name_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="string_type_info"><xsl:value-of select="string_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="tax_set_type_info"><xsl:value-of select="tax_set_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="tree_set_type_info"><xsl:value-of select="tree_set_type_info/default/constant_val"/></xsl:when>
        <xsl:when test="integer_type_info"><xsl:value-of select="integer_type_info/default/constant_val"/></xsl:when>
	</xsl:choose></xsl:variable>
	<xsl:variable name="font_style">Verdana-Plain-14</xsl:variable>
	<xsl:variable name="label_size">350,20</xsl:variable>
	<xsl:variable name="tf_size">150,20</xsl:variable>
    <xsl:variable name="cmd_append">Append to File</xsl:variable>
    <xsl:variable name="cmd_replace">Replace File</xsl:variable>

	<xsl:choose>
         <xsl:when test="./choice_type_info">
            	<hbox alignmentX="0.0" focusable="false">
					<xsl:variable name="text"><xsl:call-template name="standardize_length">
						<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
						<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
					</xsl:call-template></xsl:variable>
					<label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
            	<xsl:choose>
            		<xsl:when test="count(./choice_type_info/choices/constant_val)>4 or ./choice_type_info/choices/labile">
            			<!--combobox initclass="{$cmd_opt_combobox_class}" id="{$cmd_opt_combobox_choice}" name="{$cmd_opt_combobox_choice}" Action="{$cmd_opt_combobox_action}" ActionCommand="{$cmd_opt_combobox}" alignmentY="0.0"/-->
            			<combobox id="{$cmd_opt_combobox_choice}" name="{$cmd_opt_combobox_choice}" Action="{$cmd_opt_combobox_action}" ActionCommand="{$cmd_opt_combobox}" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
            		</xsl:when>
            		<xsl:otherwise>
            			<buttongroup>
            			<xsl:if test = "boolean(string(./choice_type_info/choices/constant_val))">
            				<xsl:for-each select = "./choice_type_info/choices/constant_val">
						<xsl:call-template name="add_radio_button">
							<xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param>
							<xsl:with-param name="number"><xsl:value-of select="position() - 1"/></xsl:with-param>
							<xsl:with-param name="font_style"><xsl:value-of select="$font_style"/></xsl:with-param>
						</xsl:call-template>					
            				</xsl:for-each>
            			</xsl:if>
            			</buttongroup>
            		</xsl:otherwise>
            	</xsl:choose>
            	</hbox>
        </xsl:when>
        <xsl:when test="(./restricted_string_type_info) or (./name_type_info)">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}" alignmentY="0.0" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_res_string}" name="{$cmd_opt_text_field_res_string}" text="{$cmd_opt_default}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}"/>
            </hbox>
        </xsl:when>
        <xsl:when test="./bool_type_info">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:choose>
					<xsl:when test = "string($cmd_opt_default)"> 
						<checkbox id="{$cmd_opt_checkbox}" name="{$cmd_opt_checkbox}" selected="{$cmd_opt_default}" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
					</xsl:when>
					<xsl:otherwise> 
						<checkbox id="{$cmd_opt_checkbox}" name="{$cmd_opt_checkbox}" selected ="false" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
					</xsl:otherwise>
				</xsl:choose>
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
            </hbox>
        </xsl:when>
        <!-- mixed command options may contain choice commands that create enums -->
        <xsl:when test="./mixed_type_info">
            <xsl:for-each select="./mixed_type_info/cmd_param">
            <xsl:sort select="@gui_order" data-type="number"/>
            <xsl:if test="position()=1">
				<xsl:call-template name="process_cmd_param">
					<xsl:with-param name="mixed">true</xsl:with-param>
					<xsl:with-param name="cmd_opt_gui_desc"><xsl:value-of select="../../gui_description"/></xsl:with-param>
					<xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param>
				</xsl:call-template>
            </xsl:if>
            <xsl:if test="position()!=1">
				<xsl:call-template name="process_cmd_param">
					<xsl:with-param name="mixed">true</xsl:with-param>
					<xsl:with-param name="cmd_opt_gui_desc"></xsl:with-param>
					<xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param>
				</xsl:call-template>
            </xsl:if>
            </xsl:for-each>
        </xsl:when>
        <xsl:when test="./outfile_type_info or ./infile_type_info" >
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_file}" name="{$cmd_opt_text_field_file}" text="{$cmd_opt_default}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}"/>
				<button text="Choose" id="choose_{$cmd_opt_choose}" Action="chooseAction_{$cmd_opt_choose}" ActionCommand = "chooseAction_{$cmd_opt_choose}" alignmentY="0.0"/>
            </hbox>
        </xsl:when>
        <xsl:when test="./output_type_info" >
			<xsl:variable name="rb_name"><xsl:text>rb_Out_cmd</xsl:text><xsl:value-of select="$cmd_name"/><xsl:text>_</xsl:text><xsl:value-of select="./phycas_impl/manipulated_var"/></xsl:variable>
			<xsl:variable name="rb_append"><xsl:value-of select="$rb_name"/><xsl:text>_Append</xsl:text></xsl:variable>
			<xsl:variable name="rb_replace"><xsl:value-of select="$rb_name"/><xsl:text>_Replace</xsl:text></xsl:variable>
			<xsl:variable name="select_file">
				<xsl:choose>
					<xsl:when test="string(./output_type_info/default/file/@path) and ./output_type_info/default/@suppress = 'false' ">true</xsl:when>
					<xsl:otherwise>false</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			<xsl:variable name="select_out">
				<xsl:choose>
					<xsl:when test="./output_type_info/default/redirect = 'Output' and ./output_type_info/default/@suppress = 'false' ">true</xsl:when>
					<xsl:otherwise>false</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			<xsl:variable name="select_plot">
				<xsl:choose>
					<xsl:when test="./output_type_info/default/redirect = 'Plot' and ./output_type_info/default/@suppress = 'false' ">true</xsl:when>
					<xsl:otherwise>false</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			<xsl:variable name="select_append">
				<xsl:choose>
					<xsl:when test="string(./output_type_info/default/file/@append)"><xsl:value-of select="./output_type_info/default/file/@append"/></xsl:when>
					<xsl:otherwise>false</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			<xsl:variable name="select_replace">
				<xsl:choose>
					<xsl:when test="string(./output_type_info/default/file/@replace)"><xsl:value-of select="./output_type_info/default/file/@replace"/></xsl:when>
					<xsl:otherwise>false</xsl:otherwise>
				</xsl:choose>
			</xsl:variable>
			<xsl:variable name="filename">
				<xsl:choose>
					<xsl:when test="string(./output_type_info/default/file/@path)"><xsl:value-of select="./output_type_info/default/file/@path"/></xsl:when>
					<xsl:otherwise></xsl:otherwise>
				</xsl:choose>
			</xsl:variable>

			<hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
				<label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="300,20" preferredSize="300,20" focusable="false"/>
				<checkbox id="{$cmd_opt_checkbox_out}_File" name="{$cmd_opt_checkbox_out}_File" selected="{$select_file}" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
				<label id="{$cmd_opt_label_out}_File" font="{$font_style}" foreground="black" text="File" alignmentY="0.0" focusable="false"/>
				<checkbox id="{$cmd_opt_checkbox_out}_Out" name="{$cmd_opt_checkbox_out}_Out" selected="{$select_out}" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
				<label id="{$cmd_opt_label_out}_Out" font="{$font_style}" foreground="black" text="Output" alignmentY="0.0"  focusable="false"/>
				<xsl:if test="./output_type_info/default/@plottable = 'true'">
					<checkbox id="{$cmd_opt_checkbox_out}_Plot" name="{$cmd_opt_checkbox_out}_Plot" selected="{$select_plot}" alignmentY="0.0" toolTipText = "{$cmd_opt_desc}"/>
					<label id="{$cmd_opt_label_out}_Plot" font="{$font_style}" foreground="black" text="Plot" alignmentY="0.0"  focusable="false"/>
				</xsl:if>
			</hbox>
            <hbox alignmentX="0.0" focusable="false">
				<label id="{$cmd_opt_label}_empty" font="{$font_style}" foreground="black" text="" alignmentY="0.0" maximumSize="75,20" preferredSize="75,20" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_out}" name="{$cmd_opt_text_field_out}" text="{$filename}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}" enabled="{$select_file}"/>
				<button text="Choose" id="choose_{$cmd_opt_choose}" name="choose_{$cmd_opt_choose}" Action="chooseAction_{$cmd_opt_choose}" ActionCommand = "chooseAction_{$cmd_opt_choose}" alignmentY="0.0" enabled="{$select_file}"/>
				<buttongroup enabled="{$select_file}">
					<radiobutton id="{$rb_append}" font="{$font_style}" name="{$rb_append}" label="Append to File" action="action_{$rb_append}" actionCommand="{$rb_append}" alignmentY="0.0" selected = "{$select_append}" enabled="{$select_file}"/>
					<radiobutton id="{$rb_replace}" font="{$font_style}" name="{$rb_replace}" label="Replace File" action="action_{$rb_replace}" actionCommand="{$rb_replace}" alignmentY="0.0" selected = "{$select_replace}" enabled="{$select_file}"/>
				</buttongroup>
            </hbox>
        </xsl:when>
        <xsl:when test="./tax_set_type_info or ./char_set_type_info or ./tree_set_type_info">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label  text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}" alignmentY="0.0" focusable="false"/>
				<scrollpane focusable="false" preferredSize="150,100" maximumSize="150,100" ><list visibleRowCount="5" id="{$cmd_opt_list_set}" name="{$cmd_opt_list_set}" Action="{$cmd_opt_list_action}" ActionCommand="{$cmd_opt_list}" toolTipText = "{$cmd_opt_desc}" /></scrollpane>
			</hbox>
			<hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="$max_gui_desc_length"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}2" alignmentY="0.0" focusable="false" maximumSize="{$label_size}" preferredSize="{$label_size}" />
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}3" alignmentY="0.0" focusable="false"/>
				<combobox id="{$cmd_opt_combobox_set}" name="{$cmd_opt_combobox_set}" Action="{$cmd_opt_combobox_action}" ActionCommand="{$cmd_opt_combobox_set}" toolTipText = "Choose a set." maximumSize="{$tf_size}" preferredSize="{$tf_size}"/>
            </hbox>
			<hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="$max_gui_desc_length"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}4" alignmentY="0.0" focusable="false" maximumSize="{$label_size}" preferredSize="{$label_size}" />
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}5" alignmentY="0.0" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_set}" name="{$cmd_opt_text_field_set}" text="" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "Enter a set."/>
            </hbox>
        </xsl:when>
        <xsl:when test="./string_type_info">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}" alignmentY="0.0" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_string}" name="{$cmd_opt_text_field_string}" text="{$cmd_opt_default}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}"/>
            </hbox>
        </xsl:when>
        <xsl:when test="./distribution_type_info">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}" alignmentY="0.0" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_string}" name="{$cmd_opt_text_field_string}" text="{$cmd_opt_default}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}"/>
				<button text="Choose" id="choose_{$cmd_opt_choose}" Action="chooseAction_{$cmd_opt_choose}" ActionCommand = "chooseAction_{$cmd_opt_choose}" alignmentY="0.0"/>
            </hbox>
        </xsl:when>
        <xsl:when test="./integer_type_info or ./double_type_info">
            <hbox alignmentX="0.0" focusable="false">
				<xsl:variable name="text"><xsl:call-template name="standardize_length">
					<xsl:with-param name="desired_length"><xsl:value-of select="$max_gui_desc_length"/></xsl:with-param>
					<xsl:with-param name="text"><xsl:value-of select="$cmd_opt_gui_desc"/></xsl:with-param>
				</xsl:call-template></xsl:variable>
                <label id="{$cmd_opt_label}" font="{$font_style}" foreground="black" text="{$text}" alignmentY="0.0" maximumSize="{$label_size}" preferredSize="{$label_size}" focusable="false"/>
				<xsl:variable name="textSpace"><xsl:call-template name="append-pad">
					<xsl:with-param name="padChar" select="' '"/>
					<xsl:with-param name="length" select="5"/>
				</xsl:call-template></xsl:variable>
				<label text="{$textSpace}" foreground="black" font="{$font_style}" id="{$cmd_opt_label_space}" alignmentY="0.0" focusable="false"/>
                <textfield id="{$cmd_opt_text_field_number}" name="{$cmd_opt_text_field_number}" text="{$cmd_opt_default}" alignmentY="0.0" maximumSize="{$tf_size}" preferredSize="{$tf_size}" toolTipText = "{$cmd_opt_desc}"/>
            </hbox>
        </xsl:when>
	</xsl:choose>
</xsl:template>

<xsl:template name="add_radio_button">
	<xsl:param name="cmd_name"></xsl:param>
     <xsl:param name="number"></xsl:param>
     <xsl:param name="font_style"></xsl:param>
	<xsl:variable name="cmd_opt_radio_button">
		<xsl:call-template name="name_radio_button"><xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param><xsl:with-param name="number"><xsl:value-of select="$number"/></xsl:with-param></xsl:call-template>
	</xsl:variable>
    <xsl:variable name="cmd_opt_radio_button_action">action_<xsl:call-template name="name_radio_button"><xsl:with-param name="cmd_name"><xsl:value-of select="$cmd_name"/></xsl:with-param></xsl:call-template>
    </xsl:variable>
    <radiobutton id="{$cmd_opt_radio_button}" font="{$font_style}" name="{$cmd_opt_radio_button}" label="{.}" action="{$cmd_opt_radio_button_action}" actionCommand="{$cmd_opt_radio_button}" alignmentY="0.0"/>
</xsl:template>

<xsl:template name="name_radio_button">
     <xsl:param name="cmd_name"></xsl:param>
     <xsl:param name="number"></xsl:param>
	<xsl:text>rb_Choice_cmd</xsl:text><xsl:value-of select="$cmd_name"/><xsl:text>_</xsl:text><xsl:value-of select="../../../phycas_impl/manipulated_var"/>
	<xsl:text>_</xsl:text>
	<xsl:value-of select ="$number"/>
</xsl:template>

<xsl:template name="add_cmd_name">
     <xsl:param name="cmd_name"></xsl:param>
	<xsl:text>cmd</xsl:text><xsl:value-of select="$cmd_name"/><xsl:text>_</xsl:text><xsl:value-of select="phycas_impl/manipulated_var"/>
</xsl:template>

<xsl:template name="calc_gui_desc_length">
	<xsl:param name="cmd_opt_gui_desc"></xsl:param>
	<xsl:param name="longest_gui_desc">0</xsl:param>
	<xsl:variable name="gui_desc_length"><xsl:value-of select="string-length($cmd_opt_gui_desc)"></xsl:value-of></xsl:variable>
	<xsl:choose>
	<xsl:when test="$gui_desc_length>$longest_gui_desc">
		<xsl:text>longer</xsl:text>
		<xsl:value-of select="string-length($cmd_opt_gui_desc)"></xsl:value-of>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>shorter</xsl:text>
		<xsl:value-of select="$longest_gui_desc"></xsl:value-of>
	</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="standardize_length">
     <xsl:param name="desired_length"></xsl:param>
     <xsl:param name="text"></xsl:param>
     <xsl:choose>
		<xsl:when test="$desired_length>string-length($text)">
			<xsl:variable name="new_text"><xsl:call-template name="append-pad">
				<xsl:with-param name="padChar" select="' '"/>
				<xsl:with-param name="padVar" select="$text"/>
				<xsl:with-param name="length" select="$desired_length"/>
			</xsl:call-template></xsl:variable>
			<xsl:value-of select="$new_text"/>
		</xsl:when>
		<xsl:otherwise>
			<xsl:value-of select="$text"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="append-pad">
  <!-- recursive template to left justify and append  -->
  <!-- the value with whatever padChar is passed in   -->
    <xsl:param name="padChar"> </xsl:param>
    <xsl:param name="padVar"/>
    <xsl:param name="length"/>
    <xsl:choose>
      <xsl:when test="string-length($padVar) &lt; $length">
        <xsl:call-template name="append-pad">
          <xsl:with-param name="padChar" select="$padChar"/>
          <xsl:with-param name="padVar" select="concat($padVar,$padChar)"/>
          <xsl:with-param name="length" select="$length"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="substring($padVar,1,$length)"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>


</xsl:stylesheet>
