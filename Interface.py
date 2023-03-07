#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st

password = "streamlit123"

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["password"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # don't store password
        else:
            st.session_state["password_correct"] = False

    if "password_correct" not in st.session_state:
        # First run, show input for password.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        return False
    elif not st.session_state["password_correct"]:
        # Password not correct, show input + error.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        st.error("😕 Password incorrect")
        return False
    else:
        # Password correct.
        return True

if check_password():

    st.set_page_config(
        page_title="Interface",
        page_icon="💻",
        )

    st.write("# Welcome to the future! 👋")

    st.sidebar.success("Select the application:")

    st.markdown(
        """
        This initiative aims to migrate applications used by the POETS INEOS engineering team to a more friendly, versatile and agile platform.
    Everyone's contribution is essential. It is a beta version, thus, if you find any abnormality, please contact the developer.
    """
    )


    st.caption('Application developed by Adilton Lopes da Silva (INEOS Polymers Engineering & Technology Support)')

else:
    


# In[ ]:




