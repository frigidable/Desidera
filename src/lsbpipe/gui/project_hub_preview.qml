import QtQuick 2.15
import QtQuick.Controls 2.15

// Isolated Observatory Studio comparison, no simulated science results.
// All project names, statuses and commands originate in the live Widgets hub.
Rectangle {
    id: root
    objectName: "observatory_qml_hub_root"
    width: 1180
    height: 740
    property bool darkMode: hubBridge.darkTheme
    property color backgroundColor: darkMode ? "#0C1424" : "#F2F5FA"
    property color surfaceColor: darkMode ? "#172238" : "#FFFFFF"
    property color raisedColor: darkMode ? "#1D2C46" : "#E9F0FA"
    property color textColor: darkMode ? "#E9F1FD" : "#1A2B45"
    property color mutedColor: darkMode ? "#AABAD1" : "#4D617A"
    property color borderColor: darkMode ? "#344862" : "#CED9E8"
    property color accentColor: darkMode ? "#65B8FA" : "#215CAA"
    property color selectionColor: darkMode ? "#1A3552" : "#E3F0FA"
    property color warnColor: darkMode ? "#EEC47A" : "#865313"
    property color warnBackground: darkMode ? "#33291D" : "#FFF3DE"
    color: backgroundColor

    Rectangle {
        id: topbar
        x: 22; y: 18; width: root.width - 44; height: 64
        radius: 12; color: root.surfaceColor
        border.width: 1; border.color: root.borderColor
        Text {
            x: 20; y: 15
            text: "◉  DESIDERA"
            font.pixelSize: 21
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 224; y: 22
            text: "OBSERVATORY STUDIO    /    PROJECT HUB"
            font.pixelSize: 12
            font.letterSpacing: 0.8
            color: root.mutedColor
        }
        Rectangle {
            width: 154; height: 27; radius: 7
            anchors.right: parent.right
            anchors.rightMargin: 20
            anchors.verticalCenter: parent.verticalCenter
            color: root.raisedColor
            Text {
                anchors.centerIn: parent
                text: "QML comparison · beta"
                font.pixelSize: 11
                color: root.mutedColor
            }
        }
    }

    Rectangle {
        id: library
        x: 22; y: 98
        width: Math.max(338, Math.min(440, root.width * 0.36))
        height: root.height - y - 56
        radius: 12
        color: root.surfaceColor
        border.width: 1
        border.color: root.borderColor

        Text {
            id: libraryHeading
            x: 20; y: 19
            text: "Galaxy projects"
            font.pixelSize: 21
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 20; y: 54
            text: hubBridge.visibleCount === hubBridge.projectCount
                ? hubBridge.projectCount + " saved galaxy projects"
                : hubBridge.visibleCount + " of " + hubBridge.projectCount + " galaxies"
            font.pixelSize: 12
            color: root.mutedColor
        }
        TextField {
            id: queryField
            objectName: "qml_hub_filter"
            x: 20; y: 87; width: library.width - 40; height: 40
            placeholderText: "Search galaxies or paths..."
            font.pixelSize: 13
            color: root.textColor
            placeholderTextColor: root.mutedColor
            selectByMouse: true
            onTextChanged: hubBridge.setFilter(text)
            background: Rectangle {
                radius: 8
                color: root.raisedColor
                border.width: queryField.activeFocus ? 2 : 1
                border.color: queryField.activeFocus ? root.accentColor : root.borderColor
            }
        }
        Text {
            x: 20; y: 144
            text: "RECENT GALAXIES"
            color: root.accentColor
            font.pixelSize: 11
            font.weight: Font.DemiBold
            font.letterSpacing: 1
        }
        ListView {
            id: galaxyList
            objectName: "qml_hub_list"
            x: 12; y: 169
            width: library.width - 24
            height: library.height - 270
            clip: true
            spacing: 8
            focus: true
            model: projectModel
            keyNavigationEnabled: true
            onCurrentItemChanged: {
                if (currentItem && activeFocus)
                    hubBridge.choose(currentItem.projectId)
            }
            Keys.onReturnPressed: hubBridge.openSelected()
            Keys.onEnterPressed: hubBridge.openSelected()
            delegate: Rectangle {
                required property int index
                required property string galaxyName
                required property string statusText
                required property string coordinates
                required property string projectId
                required property bool isAvailable
                width: galaxyList.width
                height: 93
                radius: 9
                color: hubBridge.selectedId === projectId
                    ? root.selectionColor : root.raisedColor
                border.width: hubBridge.selectedId === projectId ? 2 : 1
                border.color: hubBridge.selectedId === projectId
                    ? root.accentColor : root.borderColor
                Rectangle {
                    x: 13; y: 16; width: 40; height: 40; radius: 9
                    color: root.darkMode ? "#253955" : "#D9E8F6"
                    Text {
                        anchors.centerIn: parent
                        text: "◉"
                        font.pixelSize: 25
                        color: root.accentColor
                    }
                }
                Text {
                    x: 64; y: 13
                    width: parent.width - 79
                    text: galaxyName
                    elide: Text.ElideRight
                    color: root.textColor
                    font.pixelSize: 15
                    font.weight: Font.DemiBold
                }
                Text {
                    x: 64; y: 39
                    width: parent.width - 79
                    text: statusText
                    elide: Text.ElideRight
                    color: isAvailable ? root.accentColor : root.warnColor
                    font.pixelSize: 12
                }
                Text {
                    x: 64; y: 63
                    width: parent.width - 79
                    text: coordinates
                    elide: Text.ElideRight
                    color: root.mutedColor
                    font.pixelSize: 11
                }
                MouseArea {
                    anchors.fill: parent
                    onClicked: hubBridge.choose(projectId)
                    onDoubleClicked: {
                        hubBridge.choose(projectId)
                        hubBridge.openSelected()
                    }
                }
            }
            ScrollBar.vertical: ScrollBar { policy: ScrollBar.AsNeeded }
        }
        Text {
            anchors.centerIn: galaxyList
            width: galaxyList.width - 44
            visible: hubBridge.visibleCount === 0
            text: hubBridge.projectCount === 0
                ? "No saved galaxies yet.\nCreate a new project or import a folder."
                : "No matching galaxies.\nTry another search."
            wrapMode: Text.WordWrap
            horizontalAlignment: Text.AlignHCenter
            lineHeight: 1.4
            color: root.mutedColor
            font.pixelSize: 14
        }
        Button {
            id: newProjectButton
            objectName: "qml_hub_new_project"
            x: 20; y: library.height - 85
            width: (library.width - 48) * 0.53
            height: 45
            text: "+  New galaxy"
            onClicked: hubBridge.createProject()
            background: Rectangle {
                radius: 8
                color: newProjectButton.down ? "#16538E" : root.accentColor
            }
            contentItem: Text {
                text: newProjectButton.text
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
                font.pixelSize: 13
                font.weight: Font.DemiBold
                color: root.darkMode ? "#071629" : "#FFFFFF"
            }
        }
        Button {
            id: importButton
            objectName: "qml_hub_import_folder"
            x: newProjectButton.x + newProjectButton.width + 8
            y: newProjectButton.y
            width: library.width - x - 20
            height: 45
            text: "Import folder"
            onClicked: hubBridge.importFolder()
            background: Rectangle {
                radius: 8
                color: root.raisedColor
                border.width: 1
                border.color: root.borderColor
            }
            contentItem: Text {
                text: importButton.text
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
                font.pixelSize: 13
                color: root.textColor
            }
        }
    }

    Rectangle {
        id: inspector
        x: library.x + library.width + 16
        y: library.y
        width: root.width - x - 22
        height: library.height
        radius: 12
        color: root.surfaceColor
        border.width: 1
        border.color: root.borderColor

        Rectangle {
            x: 20; y: 20; width: 76; height: 76; radius: 11
            color: root.raisedColor
            border.width: 1
            border.color: root.borderColor
            Text {
                anchors.centerIn: parent
                text: "◉"
                font.pixelSize: 37
                color: root.accentColor
            }
        }
        Text {
            x: 114; y: 20
            width: inspector.width - 260
            text: hubBridge.selectedName
            elide: Text.ElideRight
            font.pixelSize: 25
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 114; y: 61
            width: inspector.width - 146
            text: hubBridge.selectedCoords
            wrapMode: Text.WordWrap
            font.pixelSize: 12
            color: root.mutedColor
        }
        Rectangle {
            x: inspector.width - 144; y: 21
            width: 123; height: 27; radius: 7
            color: root.warnBackground
            Text {
                anchors.centerIn: parent
                text: "SCIENCE HOLD"
                color: root.warnColor
                font.pixelSize: 12
                font.weight: Font.DemiBold
            }
        }
        Rectangle {
            id: separator
            x: 20; y: 117
            width: inspector.width - 40; height: 1
            color: root.borderColor
        }

        Text {
            x: 22; y: 136
            text: "PROJECT EVIDENCE"
            color: root.accentColor
            font.pixelSize: 11
            font.weight: Font.DemiBold
            font.letterSpacing: 1.0
        }
        Text {
            x: 22; y: 167
            text: "Preparation status"
            color: root.mutedColor
            font.pixelSize: 12
        }
        Text {
            x: 206; y: 167
            width: inspector.width - 228
            text: hubBridge.selectedStatus
            color: root.textColor
            font.pixelSize: 13
            font.weight: Font.DemiBold
            wrapMode: Text.WordWrap
        }
        Text {
            x: 22; y: 216
            text: "Survey dataset"
            color: root.mutedColor
            font.pixelSize: 12
        }
        Text {
            x: 206; y: 216
            width: inspector.width - 228
            text: hubBridge.datasetLabel
            color: root.textColor
            font.pixelSize: 12
            elide: Text.ElideRight
        }
        Text {
            x: 22; y: 265
            text: "Original FITS"
            color: root.mutedColor
            font.pixelSize: 12
        }
        Text {
            x: 206; y: 265
            width: inspector.width - 228
            text: hubBridge.fitsLabel
            color: root.textColor
            font.pixelSize: 12
            wrapMode: Text.WordWrap
        }
        Rectangle {
            x: 20; y: 328; width: inspector.width - 40
            height: 1; color: root.borderColor
        }
        Rectangle {
            id: nextBox
            objectName: "qml_hub_next_action_card"
            x: 20; y: 351; width: inspector.width - 40
            height: 147; radius: 10
            color: root.raisedColor
            border.width: 1; border.color: root.borderColor
            Text {
                x: 16; y: 15
                text: "RECOMMENDED NEXT STEP"
                color: root.accentColor
                font.pixelSize: 11
                font.weight: Font.DemiBold
                font.letterSpacing: 0.7
            }
            Text {
                x: 16; y: 43
                width: nextBox.width - 32
                text: hubBridge.nextAction
                color: root.textColor
                font.pixelSize: 15
                font.weight: Font.DemiBold
                wrapMode: Text.WordWrap
            }
            Text {
                x: 16; y: 104
                width: nextBox.width - 32
                text: "Based on saved project state · no science adopted"
                color: root.mutedColor
                font.pixelSize: 11
                wrapMode: Text.WordWrap
            }
        }
        Button {
            id: openButton
            objectName: "qml_hub_open_project"
            x: 20; y: inspector.height - 81
            width: inspector.width - 40; height: 47
            text: hubBridge.selectedId.length === 0 ? "Select a galaxy"
                : hubBridge.selectedStatus.indexOf("Missing") >= 0
                  ? "Relink original project..." : "Open galaxy workspace"
            enabled: hubBridge.selectedId.length !== 0
            onClicked: hubBridge.openSelected()
            background: Rectangle {
                radius: 8
                color: openButton.enabled ? root.accentColor : root.raisedColor
            }
            contentItem: Text {
                text: openButton.text
                color: openButton.enabled
                    ? (root.darkMode ? "#071629" : "#FFFFFF") : root.mutedColor
                font.pixelSize: 14
                font.weight: Font.DemiBold
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
            }
        }
    }

    Text {
        x: 28; y: root.height - 34
        text: "Original FITS are read-only · catalogue associations are not scientific adoption"
        color: root.mutedColor
        font.pixelSize: 11
    }
    Text {
        anchors.right: parent.right
        anchors.rightMargin: 28
        y: root.height - 34
        text: "Qt Quick experimental preview"
        color: root.mutedColor
        font.pixelSize: 11
    }
}
