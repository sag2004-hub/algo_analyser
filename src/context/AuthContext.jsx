// src/context/AuthContext.jsx
import React, { createContext, useContext, useEffect, useState } from 'react';
import { auth, db, serverTs } from '../firebase';
import { onAuthStateChanged, signOut } from 'firebase/auth';
import { doc, setDoc } from 'firebase/firestore';

// Create context and hook
const AuthContext = createContext(null);
export const useAuth = () => useContext(AuthContext);

// Function to upsert user data in Firestore
async function upsertUserDoc(user) {
  if (!user) return;
  const { uid, email, displayName, photoURL, providerData } = user;
  const providerId = providerData?.[0]?.providerId ?? 'password';

  await setDoc(
    doc(db, 'users', uid),
    {
      uid,
      email: email ?? '',
      displayName: displayName ?? '',
      photoURL: photoURL ?? '',
      providerId,
      lastLoginAt: serverTs(),
    },
    { merge: true }
  );
}

// AuthProvider component
export function AuthProvider({ children }) {
  const [user, setUser] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const unsub = onAuthStateChanged(auth, async (u) => {
      setUser(u);
      if (u) await upsertUserDoc(u);
      setLoading(false);
    });

    return () => unsub();
  }, []);

  const value = {
    user,
    loading,
    signOutUser: () => signOut(auth),
  };

  // Prevent rendering children until auth state is loaded
  return (
    <AuthContext.Provider value={value}>
      {!loading && children}
    </AuthContext.Provider>
  );
}
