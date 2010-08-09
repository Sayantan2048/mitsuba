/*
* Copyright 2006 Sony Computer Entertainment Inc.
*
* Licensed under the MIT Open Source License, for details please see license.txt or the website
* http://www.opensource.org/licenses/mit-license.php
*
*/ 

#ifndef __domLibrary_physics_scenes_h__
#define __domLibrary_physics_scenes_h__

#include <dae/daeDocument.h>
#include <dom/domTypes.h>
#include <dom/domElements.h>

#include <dom/domAsset.h>
#include <dom/domPhysics_scene.h>
#include <dom/domExtra.h>
class DAE;

/**
 * The library_physics_scenes element declares a module of physics_scene elements.
 */
class domLibrary_physics_scenes : public daeElement
{
public:
	virtual COLLADA_TYPE::TypeEnum getElementType() const { return COLLADA_TYPE::LIBRARY_PHYSICS_SCENES; }
	static daeInt ID() { return 725; }
	virtual daeInt typeID() const { return ID(); }
protected:  // Attributes
/**
 *  The id attribute is a text string containing the unique identifier of
 * this element.  This value must be unique within the instance document.
 * Optional attribute. 
 */
	xsID attrId;
/**
 *  The name attribute is the text string name of this element. Optional attribute.
 */
	xsNCName attrName;

protected:  // Elements
/**
 *  The library_physics_scenes element may contain an asset element.  @see
 * domAsset
 */
	domAssetRef elemAsset;
/**
 *  There must be at least one physics_scene element.  @see domPhysics_scene
 */
	domPhysics_scene_Array elemPhysics_scene_array;
/**
 *  The extra element may appear any number of times.  @see domExtra
 */
	domExtra_Array elemExtra_array;

public:	//Accessors and Mutators
	/**
	 * Gets the id attribute.
	 * @return Returns a xsID of the id attribute.
	 */
	xsID getId() const { return attrId; }
	/**
	 * Sets the id attribute.
	 * @param atId The new value for the id attribute.
	 */
	void setId( xsID atId ) { *(daeStringRef*)&attrId = atId; _validAttributeArray[0] = true; 
		if( _document != NULL ) _document->changeElementID( this, attrId );
	}

	/**
	 * Gets the name attribute.
	 * @return Returns a xsNCName of the name attribute.
	 */
	xsNCName getName() const { return attrName; }
	/**
	 * Sets the name attribute.
	 * @param atName The new value for the name attribute.
	 */
	void setName( xsNCName atName ) { *(daeStringRef*)&attrName = atName; _validAttributeArray[1] = true; }

	/**
	 * Gets the asset element.
	 * @return a daeSmartRef to the asset element.
	 */
	const domAssetRef getAsset() const { return elemAsset; }
	/**
	 * Gets the physics_scene element array.
	 * @return Returns a reference to the array of physics_scene elements.
	 */
	domPhysics_scene_Array &getPhysics_scene_array() { return elemPhysics_scene_array; }
	/**
	 * Gets the physics_scene element array.
	 * @return Returns a constant reference to the array of physics_scene elements.
	 */
	const domPhysics_scene_Array &getPhysics_scene_array() const { return elemPhysics_scene_array; }
	/**
	 * Gets the extra element array.
	 * @return Returns a reference to the array of extra elements.
	 */
	domExtra_Array &getExtra_array() { return elemExtra_array; }
	/**
	 * Gets the extra element array.
	 * @return Returns a constant reference to the array of extra elements.
	 */
	const domExtra_Array &getExtra_array() const { return elemExtra_array; }
protected:
	/**
	 * Constructor
	 */
	domLibrary_physics_scenes(DAE& dae) : daeElement(dae), attrId(), attrName(), elemAsset(), elemPhysics_scene_array(), elemExtra_array() {}
	/**
	 * Destructor
	 */
	virtual ~domLibrary_physics_scenes() {}
	/**
	 * Overloaded assignment operator
	 */
	virtual domLibrary_physics_scenes &operator=( const domLibrary_physics_scenes &cpy ) { (void)cpy; return *this; }

public: // STATIC METHODS
	/**
	 * Creates an instance of this class and returns a daeElementRef referencing it.
	 * @return a daeElementRef referencing an instance of this object.
	 */
	static DLLSPEC daeElementRef create(DAE& dae);
	/**
	 * Creates a daeMetaElement object that describes this element in the meta object reflection framework.
	 * If a daeMetaElement already exists it will return that instead of creating a new one. 
	 * @return A daeMetaElement describing this COLLADA element.
	 */
	static DLLSPEC daeMetaElement* registerElement(DAE& dae);
};


#endif
