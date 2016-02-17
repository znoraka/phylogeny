#ifndef HEAP_H
#define HEAP_H

#include <iostream>
#include <vector>

template<typename T>
class Heap
{
public:
  Heap(int s, std::function<bool(T a, T b)> comp = [&](T a, T b){return a < b;});
  void add(T v);
  void pop();
  T top();
  void display();
  std::vector<T> sort();
  int size() const;

private:
  int s;
  // T* heap;
  std::vector<T> heap;
  std::function<bool(T a, T b)> comp;

  void up();
  int father(int i);
  int leftSon(int i);
  int rightSon(int i);
  bool fatherGreater(int i);
  void switchWithFather(int i);
  void down();
  int smallerSon(int i);
};

template<typename T>
Heap<T>::Heap(int size, std::function<bool(T a, T b)> comp)
{
  // heap = new T[size];
  heap.resize(size + 1);
  this->comp = comp;
  this->s = 1;
}

template<typename T>
int Heap<T>::father(int i)
{
  if(i/2 < 1) {
    return 1;
  }
  return i/2;
}

template<typename T>
int Heap<T>::leftSon(int i)
{
  return i*2;
}

template<typename T>
int Heap<T>::rightSon(int i)
{
  return i*2+1;
}

template<typename T>
bool Heap<T>::fatherGreater(int i)
{
  return comp(heap[i], heap[father(i)]);
  // return *heap[i] < *heap[father(i)];
}

template<typename T>
void Heap<T>::switchWithFather(int i)
{
  T temp = heap[i];
  heap[i] = heap[father(i)];
  heap[father(i)] = temp;
}

template<typename T>
void Heap<T>::up()
{
  for (int i = s; i > 0; i /= 2)
    {
      if(fatherGreater(i))
        {
	  switchWithFather(i);
        }
      else
        {
	  return;
        }
    }
}

template<typename T>
void Heap<T>::add(T v)
{
  heap[s] = v;
  up();
  s++;
}

template<typename T>
int Heap<T>::smallerSon(int i)
{
  if(rightSon(i) > s)
    {
      return leftSon(i);
    }

  if(comp(heap[leftSon(i)], heap[rightSon(i)]))
  // if(*heap[leftSon(i)] < *heap[rightSon(i)])
    {
      return leftSon(i);
    }
  else
    {
      return rightSon(i);
    }
}

template<typename T>
void Heap<T>::down()
{
  bool b = true;
  int i = 1;

  while(b)
    {
      b = false;

      int ss = smallerSon(i);

      if(ss > s)
	return;

      if(comp(heap[ss], heap[i]))
      // if(*heap[ss] < *heap[i])
        {
	  switchWithFather(ss);
	  b = true;
	  i = ss;
        }
    }
}

template<typename T>
void Heap<T>::pop()
{
  s--;
  heap[1] = heap[s];
  down();
}

template<typename T>
T Heap<T>::top()
{
  return heap[1];
}

template<typename T>
std::vector<T> Heap<T>::sort()
{
  // T* sorted = new T[s];
  std::vector<T> sorted;
  int n = s;

  for (int i = 0; i < n-1; ++i)
    {
      // sorted[i] = top();
      sorted.push_back(top());
      pop();

    }
  return sorted;
}

template<typename T>
int Heap<T>::size() const
{
  return s - 1;
}

template<typename T>
void Heap<T>::display()
{
  for (int i = 1; i < s; i++)
    {
      std::cout << heap[i]->value << "  ";
    }
  std::cout << std::endl;
}

#endif // HEAP_H
